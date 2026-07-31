use noodles::{
    core::Position,
    vcf::{
        self,
        header::record::value::{
            Map,
            map::{
                AlternativeAllele, Contig, Format, Info,
                format::{Number as FormatNumber, Type as FormatType},
                info::{Number as InfoNumber, Type as InfoType},
            },
        },
        variant::record_buf::{
            AlternateBases, Info as RecordInfo, Samples,
            info::field::{Value as InfoValue, value::Array as InfoArray},
            samples::{
                Keys,
                sample::{Value as SampleValue, value::Array as SampleArray},
            },
        },
    },
};

use crate::{Allele, CallError, LikelihoodSite, Ploidy, Result, SampleEvidence, SampleLikelihood};

const QUALITY_SUM: &str = "QS";
const PL: &str = "PL";
const DP: &str = "DP";
const AD: &str = "AD";

#[derive(Clone, Debug)]
pub struct LikelihoodVcfSchema {
    header: vcf::Header,
}

impl LikelihoodVcfSchema {
    pub fn new(
        references: impl IntoIterator<Item = (impl AsRef<[u8]>, u64)>,
        sample_names: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<Self> {
        let mut builder = vcf::Header::builder()
            .add_alternative_allele(
                "*",
                Map::<AlternativeAllele>::new("Represents allele(s) other than observed."),
            )
            .add_info(
                QUALITY_SUM,
                Map::<Info>::new(
                    InfoNumber::ReferenceAlternateBases,
                    InfoType::Float,
                    "Auxiliary tag used for calling",
                ),
            )
            .add_format(
                PL,
                Map::<Format>::new(
                    FormatNumber::Samples,
                    FormatType::Integer,
                    "List of Phred-scaled genotype likelihoods",
                ),
            )
            .add_format(
                DP,
                Map::<Format>::new(
                    FormatNumber::Count(1),
                    FormatType::Integer,
                    "Number of high-quality bases",
                ),
            )
            .add_format(
                AD,
                Map::<Format>::new(
                    FormatNumber::ReferenceAlternateBases,
                    FormatType::Integer,
                    "Allelic depths (high-quality bases)",
                ),
            )
            .add_format(
                QUALITY_SUM,
                Map::<Format>::new(
                    FormatNumber::ReferenceAlternateBases,
                    FormatType::Integer,
                    "Phred-score allele quality sum used by calling",
                ),
            );

        let mut reference_names = Vec::new();
        for (name, length) in references {
            let name = std::str::from_utf8(name.as_ref())
                .map_err(|_| invalid("reference name is not UTF-8"))?;
            if reference_names.iter().any(|value| value == name) {
                return Err(invalid(format!("duplicate reference name: {name}")));
            }
            reference_names.push(name.to_owned());
            let mut contig = Map::<Contig>::new();
            *contig.length_mut() = Some(
                usize::try_from(length)
                    .map_err(|_| invalid("reference length exceeds VCF range"))?,
            );
            builder = builder.add_contig(name, contig);
        }
        let mut samples = Vec::new();
        for sample_name in sample_names {
            let sample_name = sample_name.as_ref();
            if samples.iter().any(|value| value == sample_name) {
                return Err(invalid(format!("duplicate sample name: {sample_name}")));
            }
            samples.push(sample_name.to_owned());
            builder = builder.add_sample_name(sample_name);
        }

        Self::from_header(builder.build())
    }

    pub fn from_header(header: vcf::Header) -> Result<Self> {
        if header.contigs().is_empty() {
            return Err(invalid("header has no contigs"));
        }
        if header.sample_names().is_empty() {
            return Err(invalid("header has no samples"));
        }
        require_info(
            &header,
            QUALITY_SUM,
            InfoNumber::ReferenceAlternateBases,
            InfoType::Float,
        )?;
        require_format(
            &header,
            PL,
            FormatNumber::Samples,
            FormatType::Integer,
            true,
        )?;
        require_format(
            &header,
            DP,
            FormatNumber::Count(1),
            FormatType::Integer,
            false,
        )?;
        require_format(
            &header,
            AD,
            FormatNumber::ReferenceAlternateBases,
            FormatType::Integer,
            false,
        )?;
        require_format(
            &header,
            QUALITY_SUM,
            FormatNumber::ReferenceAlternateBases,
            FormatType::Integer,
            false,
        )?;
        Ok(Self { header })
    }

    pub fn header(&self) -> &vcf::Header {
        &self.header
    }

    pub fn encode_likelihood(&self, site: &LikelihoodSite) -> Result<vcf::variant::RecordBuf> {
        if site.samples().len() != self.header.sample_names().len() {
            return Err(invalid("record sample count differs from the header"));
        }
        let (reference_name, _) = self
            .header
            .contigs()
            .get_index(site.reference_sequence_id())
            .ok_or_else(|| invalid("reference sequence ID is absent from the header"))?;
        let position = site
            .position()
            .checked_add(1)
            .and_then(|position| usize::try_from(position).ok())
            .and_then(|position| Position::try_from(position).ok())
            .ok_or_else(|| invalid("position exceeds VCF range"))?;
        let alternate_bases = AlternateBases::from(
            site.alternates()
                .iter()
                .map(|allele| allele_text(allele).map(str::to_owned))
                .collect::<Result<Vec<_>>>()?,
        );
        let info = RecordInfo::from_iter([(
            QUALITY_SUM.to_owned(),
            Some(InfoValue::Array(InfoArray::Float(
                site.allele_quality_sums()
                    .iter()
                    .copied()
                    .map(Some)
                    .collect(),
            ))),
        )]);
        let keys = Keys::from_iter([
            PL.to_owned(),
            DP.to_owned(),
            AD.to_owned(),
            QUALITY_SUM.to_owned(),
        ]);
        let values = site
            .samples()
            .iter()
            .map(|sample| {
                let evidence = sample.evidence();
                Ok(vec![
                    sample
                        .phred_likelihoods()
                        .map(|values| checked_array(values, PL))
                        .transpose()?,
                    Some(SampleValue::Integer(checked_integer(evidence.depth(), DP)?)),
                    Some(checked_array(evidence.allele_depths(), AD)?),
                    Some(checked_array(evidence.allele_quality_sums(), QUALITY_SUM)?),
                ])
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(reference_name)
            .set_variant_start(position)
            .set_reference_bases(allele_text(site.reference())?)
            .set_alternate_bases(alternate_bases)
            .set_info(info)
            .set_samples(Samples::new(keys, values))
            .build())
    }

    pub fn decode_likelihood(&self, record: &vcf::variant::RecordBuf) -> Result<LikelihoodSite> {
        let reference_sequence_id = self
            .header
            .contigs()
            .get_index_of(record.reference_sequence_name())
            .ok_or_else(|| invalid("record contig is absent from the header"))?;
        let position = record
            .variant_start()
            .map(usize::from)
            .and_then(|position| position.checked_sub(1))
            .and_then(|position| u64::try_from(position).ok())
            .ok_or_else(|| invalid("record position is missing or invalid"))?;
        let reference = Allele::new(record.reference_bases().as_bytes())?;
        let alternates = record
            .alternate_bases()
            .as_ref()
            .iter()
            .map(|allele| Allele::new(allele.as_bytes()))
            .collect::<Result<Vec<_>>>()?;
        if alternates.is_empty() {
            return Err(invalid("record has no alternate alleles"));
        }
        let allele_count = alternates.len() + 1;
        let allele_quality_sums = info_quality_sums(record, allele_count)?;
        if record.samples().values().count() != self.header.sample_names().len() {
            return Err(invalid("record sample count differs from the header"));
        }
        if record.samples().select(PL).is_none() {
            return Err(invalid("record has no FORMAT/PL field"));
        }

        let samples = (0..self.header.sample_names().len())
            .map(|index| decode_sample(record.samples(), index, allele_count))
            .collect::<Result<Vec<_>>>()?;
        LikelihoodSite::new(
            reference_sequence_id,
            position,
            reference,
            alternates,
            allele_quality_sums,
            samples,
        )
    }
}

fn require_info(header: &vcf::Header, key: &str, number: InfoNumber, ty: InfoType) -> Result<()> {
    let value = header
        .infos()
        .get(key)
        .ok_or_else(|| invalid(format!("header has no INFO/{key} definition")))?;
    if value.number() != number || value.ty() != ty {
        return Err(invalid(format!(
            "header has an incompatible INFO/{key} definition"
        )));
    }
    Ok(())
}

fn require_format(
    header: &vcf::Header,
    key: &str,
    number: FormatNumber,
    ty: FormatType,
    required: bool,
) -> Result<()> {
    let Some(value) = header.formats().get(key) else {
        if required {
            return Err(invalid(format!("header has no FORMAT/{key} definition")));
        }
        return Ok(());
    };
    if value.number() != number || value.ty() != ty {
        return Err(invalid(format!(
            "header has an incompatible FORMAT/{key} definition"
        )));
    }
    Ok(())
}

fn decode_sample(samples: &Samples, index: usize, allele_count: usize) -> Result<SampleLikelihood> {
    let phred_likelihoods = sample_integer_array(samples, PL, index)?;
    let diploid_count = allele_count
        .checked_add(1)
        .and_then(|value| allele_count.checked_mul(value))
        .map(|value| value / 2)
        .ok_or_else(|| invalid("record allele count exceeds the supported range"))?;
    let ploidy = match phred_likelihoods.as_deref() {
        Some(values) if values.len() == allele_count => Ploidy::new(1)?,
        Some(values) if values.len() == diploid_count => Ploidy::new(2)?,
        Some(_) => return Err(invalid("FORMAT/PL has an unsupported genotype count")),
        None => Ploidy::new(2)?,
    };
    let allele_depths =
        sample_integer_array(samples, AD, index)?.unwrap_or_else(|| vec![0; allele_count]);
    let allele_quality_sums =
        sample_integer_array(samples, QUALITY_SUM, index)?.unwrap_or_else(|| vec![0; allele_count]);
    if allele_depths.len() != allele_count || allele_quality_sums.len() != allele_count {
        return Err(invalid(
            "FORMAT allele arrays do not match the record alleles",
        ));
    }
    let depth = match sample_integer(samples, DP, index)? {
        Some(depth) => depth,
        None => allele_depths.iter().try_fold(0u32, |total, &value| {
            total
                .checked_add(value)
                .ok_or_else(|| invalid("FORMAT/AD sum exceeds the supported range"))
        })?,
    };
    let evidence = SampleEvidence::new(depth, allele_depths, allele_quality_sums)?;
    SampleLikelihood::new(
        ploidy,
        phred_likelihoods.map(Vec::into_boxed_slice),
        evidence,
    )
}

fn info_quality_sums(record: &vcf::variant::RecordBuf, allele_count: usize) -> Result<Vec<f32>> {
    let value = record
        .info()
        .get(QUALITY_SUM)
        .flatten()
        .ok_or_else(|| invalid("record has no INFO/QS value"))?;
    let InfoValue::Array(InfoArray::Float(values)) = value else {
        return Err(invalid("INFO/QS is not a float array"));
    };
    if values.len() != allele_count {
        return Err(invalid("INFO/QS does not match the record alleles"));
    }
    values
        .iter()
        .map(|value| value.ok_or_else(|| invalid("INFO/QS contains a missing value")))
        .collect()
}

fn sample_integer(samples: &Samples, key: &str, index: usize) -> Result<Option<u32>> {
    let Some(series) = samples.select(key) else {
        return Ok(None);
    };
    match series.get(index) {
        Some(None) => Ok(None),
        Some(Some(SampleValue::Integer(value))) => u32::try_from(*value)
            .map(Some)
            .map_err(|_| invalid(format!("FORMAT/{key} contains a negative integer"))),
        Some(Some(_)) => Err(invalid(format!("FORMAT/{key} is not an integer"))),
        None => Err(invalid("record sample fields are truncated")),
    }
}

fn sample_integer_array(samples: &Samples, key: &str, index: usize) -> Result<Option<Vec<u32>>> {
    let Some(series) = samples.select(key) else {
        return Ok(None);
    };
    match series.get(index) {
        Some(None) => Ok(None),
        Some(Some(SampleValue::Array(SampleArray::Integer(values)))) => values
            .iter()
            .map(|value| {
                let value = value
                    .ok_or_else(|| invalid(format!("FORMAT/{key} contains a missing value")))?;
                u32::try_from(value)
                    .map_err(|_| invalid(format!("FORMAT/{key} contains a negative integer")))
            })
            .collect::<Result<Vec<_>>>()
            .map(Some),
        Some(Some(_)) => Err(invalid(format!("FORMAT/{key} is not an integer array"))),
        None => Err(invalid("record sample fields are truncated")),
    }
}

fn checked_array(values: &[u32], key: &str) -> Result<SampleValue> {
    values
        .iter()
        .map(|&value| checked_integer(value, key).map(Some))
        .collect::<Result<Vec<_>>>()
        .map(SampleArray::Integer)
        .map(SampleValue::Array)
}

fn checked_integer(value: u32, key: &str) -> Result<i32> {
    i32::try_from(value).map_err(|_| invalid(format!("FORMAT/{key} exceeds the VCF integer range")))
}

fn allele_text(allele: &Allele) -> Result<&str> {
    std::str::from_utf8(allele.as_bytes()).map_err(|_| invalid("allele is not UTF-8"))
}

fn invalid(message: impl Into<String>) -> CallError {
    CallError::InvalidLikelihoodVariant(message.into())
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use noodles::vcf::variant::io::Write as _;

    use super::*;

    #[test]
    fn decodes_bcftools_1_24_likelihood_record() {
        let data = include_bytes!("../tests/golden/bcftools-1.24-likelihood.vcf");
        let mut reader = vcf::io::Reader::new(&data[..]);
        let header = reader.read_header().unwrap();
        let schema = LikelihoodVcfSchema::from_header(header.clone()).unwrap();
        let mut record = vcf::variant::RecordBuf::default();
        assert_ne!(reader.read_record_buf(&header, &mut record).unwrap(), 0);

        let site = schema.decode_likelihood(&record).unwrap();

        assert_eq!(site.reference_sequence_id(), 0);
        assert_eq!(site.position(), 0);
        assert_eq!(site.reference().as_bytes(), b"A");
        assert_eq!(
            site.alternates()
                .iter()
                .map(Allele::as_bytes)
                .collect::<Vec<_>>(),
            [b"G".as_slice(), b"C".as_slice(), b"<*>".as_slice()]
        );
        assert_eq!(site.allele_quality_sums(), &[1.0, 1.0, 1.0, 0.0]);
        assert_eq!(site.samples().len(), 3);
        assert_eq!(site.samples()[0].evidence().depth(), 1);
        assert_eq!(site.samples()[1].evidence().allele_depths(), &[0, 0, 1, 0]);
        assert_eq!(
            site.samples()[2].evidence().allele_quality_sums(),
            &[0, 40, 0, 0]
        );
        assert_eq!(
            site.samples()[1].phred_likelihoods(),
            Some(&[40, 40, 40, 3, 3, 0, 40, 40, 3, 40][..])
        );
        let called = crate::MultiallelicCaller::default().call(&site).unwrap();
        assert_eq!(called.allele_counts(), &[3, 2, 1]);
        assert!((called.quality().unwrap() - 15.6934).abs() < 1e-4);

        let encoded = schema.encode_likelihood(&site).unwrap();
        assert_eq!(schema.decode_likelihood(&encoded).unwrap(), site);

        let mut vcf_data = Vec::new();
        let mut vcf_writer = vcf::io::Writer::new(&mut vcf_data);
        vcf_writer.write_header(schema.header()).unwrap();
        vcf_writer
            .write_variant_record(schema.header(), &encoded)
            .unwrap();
        assert_eq!(
            std::str::from_utf8(&vcf_data)
                .unwrap()
                .lines()
                .last()
                .unwrap(),
            "chr1\t1\t.\tA\tG,C,<*>\t.\t.\tQS=1,1,1,0\tPL:DP:AD:QS\t0,3,40,3,40,40,3,40,40,40:1:1,0,0,0:40,0,0,0\t40,40,40,3,3,0,40,40,3,40:1:0,0,1,0:0,0,40,0\t40,3,0,40,3,40,40,3,40,40:1:0,1,0,0:0,40,0,0"
        );
        let mut vcf_reader = vcf::io::Reader::new(&vcf_data[..]);
        let vcf_header = vcf_reader.read_header().unwrap();
        let vcf_schema = LikelihoodVcfSchema::from_header(vcf_header.clone()).unwrap();
        let mut vcf_record = vcf::variant::RecordBuf::default();
        vcf_reader
            .read_record_buf(&vcf_header, &mut vcf_record)
            .unwrap();
        assert_eq!(vcf_schema.decode_likelihood(&vcf_record).unwrap(), site);

        let mut bcf_data = Vec::new();
        let mut bcf_writer = noodles::bcf::io::Writer::from(&mut bcf_data);
        bcf_writer.write_header(schema.header()).unwrap();
        bcf_writer
            .write_variant_record(schema.header(), &encoded)
            .unwrap();
        let mut bcf_reader = noodles::bcf::io::Reader::from(&bcf_data[..]);
        let bcf_header = bcf_reader.read_header().unwrap();
        let bcf_schema = LikelihoodVcfSchema::from_header(bcf_header.clone()).unwrap();
        let mut bcf_record = vcf::variant::RecordBuf::default();
        bcf_reader
            .read_record_buf(&bcf_header, &mut bcf_record)
            .unwrap();
        assert_eq!(bcf_schema.decode_likelihood(&bcf_record).unwrap(), site);
    }

    #[test]
    fn builds_a_checked_likelihood_header() {
        let references = [(b"chr1".as_slice(), 5)];
        let schema = LikelihoodVcfSchema::new(references, ["s1", "s2"]).unwrap();

        assert_eq!(schema.header().contigs().len(), 1);
        assert_eq!(schema.header().sample_names().len(), 2);
        assert!(schema.header().infos().contains_key(QUALITY_SUM));
        assert!(schema.header().formats().contains_key(PL));
        assert!(schema.header().formats().contains_key(DP));
        assert!(schema.header().formats().contains_key(AD));
    }

    #[test]
    fn rejects_incompatible_likelihood_schema() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let data = std::fs::read(fixtures.join("bcftools-1.24-likelihood.vcf")).unwrap();
        let mut reader = vcf::io::Reader::new(&data[..]);
        let mut header = reader.read_header().unwrap();
        header.infos_mut().shift_remove(QUALITY_SUM);

        assert_eq!(
            LikelihoodVcfSchema::from_header(header).unwrap_err(),
            invalid("header has no INFO/QS definition")
        );
    }

    #[test]
    fn preserves_mixed_ploidy_likelihoods_in_bcf() {
        let schema =
            LikelihoodVcfSchema::new([(b"chr1".as_slice(), 5)], ["haploid", "diploid"]).unwrap();
        let site = LikelihoodSite::new(
            0,
            0,
            Allele::new(&b"A"[..]).unwrap(),
            [Allele::new(&b"G"[..]).unwrap()],
            [1.0, 1.0],
            [
                SampleLikelihood::observed(
                    Ploidy::new(1).unwrap(),
                    [0, 40],
                    SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
                )
                .unwrap(),
                SampleLikelihood::observed(
                    Ploidy::new(2).unwrap(),
                    [40, 3, 0],
                    SampleEvidence::new(1, [0, 1], [0, 40]).unwrap(),
                )
                .unwrap(),
            ],
        )
        .unwrap();
        let encoded = schema.encode_likelihood(&site).unwrap();
        let mut data = Vec::new();
        let mut writer = noodles::bcf::io::Writer::from(&mut data);
        writer.write_header(schema.header()).unwrap();
        writer
            .write_variant_record(schema.header(), &encoded)
            .unwrap();
        let mut reader = noodles::bcf::io::Reader::from(&data[..]);
        let header = reader.read_header().unwrap();
        let decoded_schema = LikelihoodVcfSchema::from_header(header.clone()).unwrap();
        let mut record = vcf::variant::RecordBuf::default();
        reader.read_record_buf(&header, &mut record).unwrap();

        assert_eq!(decoded_schema.decode_likelihood(&record).unwrap(), site);
    }
}
