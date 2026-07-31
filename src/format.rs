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
        variant::record::samples::series::value::genotype::Phasing,
        variant::record_buf::{
            AlternateBases, Info as RecordInfo, Samples,
            info::field::{Value as InfoValue, value::Array as InfoArray},
            samples::{
                Keys,
                sample::{
                    Value as SampleValue,
                    value::{Array as SampleArray, Genotype, genotype::Allele as GenotypeAllele},
                },
            },
        },
    },
};

use crate::{
    Allele, CallError, CalledSite, IndelSummary, LikelihoodSite, Ploidy, Result, SampleEvidence,
    SampleLikelihood,
};

const QUALITY_SUM: &str = "QS";
const INDEL: &str = "INDEL";
const IDV: &str = "IDV";
const IMF: &str = "IMF";
const PL: &str = "PL";
const DP: &str = "DP";
const AD: &str = "AD";
const AC: &str = "AC";
const AN: &str = "AN";
const GT: &str = "GT";
const GQ: &str = "GQ";
const GP: &str = "GP";

#[derive(Clone, Debug)]
pub struct LikelihoodVcfSchema {
    header: vcf::Header,
}

#[derive(Clone, Debug)]
pub struct CalledVcfSchema {
    header: vcf::Header,
}

impl LikelihoodVcfSchema {
    pub fn new(
        references: impl IntoIterator<Item = (impl AsRef<[u8]>, u64)>,
        sample_names: impl IntoIterator<Item = impl AsRef<str>>,
    ) -> Result<Self> {
        let mut builder = vcf::Header::builder()
            .set_file_format(vcf::header::FileFormat::new(4, 2))
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
            .add_info(
                INDEL,
                Map::<Info>::new(
                    InfoNumber::Count(0),
                    InfoType::Flag,
                    "Indicates that the variant is an INDEL",
                ),
            )
            .add_info(
                IDV,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Integer,
                    "Maximum number of raw reads supporting an indel",
                ),
            )
            .add_info(
                IMF,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Maximum fraction of raw reads supporting an indel",
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
        require_optional_info(&header, INDEL, InfoNumber::Count(0), InfoType::Flag)?;
        require_optional_info(&header, IDV, InfoNumber::Count(1), InfoType::Integer)?;
        require_optional_info(&header, IMF, InfoNumber::Count(1), InfoType::Float)?;
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
        let mut info = RecordInfo::from_iter([(
            QUALITY_SUM.to_owned(),
            Some(InfoValue::Array(InfoArray::Float(
                site.allele_quality_sums()
                    .iter()
                    .copied()
                    .map(Some)
                    .collect(),
            ))),
        )]);
        insert_indel_info(&self.header, &mut info, site.indel_summary())?;
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
        let site = LikelihoodSite::new(
            reference_sequence_id,
            position,
            reference,
            alternates,
            allele_quality_sums,
            samples,
        )?;
        match decode_indel_summary(record)? {
            Some(summary) => Ok(site.with_indel_summary(summary)),
            None => Ok(site),
        }
    }
}

impl CalledVcfSchema {
    pub fn from_likelihood(schema: &LikelihoodVcfSchema) -> Self {
        Self::from_likelihood_inner(schema, true)
    }

    pub fn from_consensus_likelihood(schema: &LikelihoodVcfSchema) -> Self {
        Self::from_likelihood_inner(schema, false)
    }

    fn from_likelihood_inner(schema: &LikelihoodVcfSchema, include_probabilities: bool) -> Self {
        let mut header = schema.header.clone();
        header.infos_mut().shift_remove(QUALITY_SUM);
        header.infos_mut().insert(
            AC.to_owned(),
            Map::<Info>::new(
                InfoNumber::AlternateBases,
                InfoType::Integer,
                "Allele count in genotypes for each ALT allele, in the same order as listed",
            ),
        );
        header.infos_mut().insert(
            AN.to_owned(),
            Map::<Info>::new(
                InfoNumber::Count(1),
                InfoType::Integer,
                "Total number of alleles in called genotypes",
            ),
        );
        header.formats_mut().insert(
            GT.to_owned(),
            Map::<Format>::new(FormatNumber::Count(1), FormatType::String, "Genotype"),
        );
        header.formats_mut().insert(
            GQ.to_owned(),
            Map::<Format>::new(
                FormatNumber::Count(1),
                FormatType::Integer,
                "Phred-scaled Genotype Quality",
            ),
        );
        if include_probabilities {
            header.formats_mut().insert(
                GP.to_owned(),
                Map::<Format>::new(
                    FormatNumber::Samples,
                    FormatType::Float,
                    "Genotype posterior probabilities in the range 0 to 1",
                ),
            );
        } else {
            header.formats_mut().shift_remove(GP);
        }
        Self { header }
    }

    pub fn header(&self) -> &vcf::Header {
        &self.header
    }

    pub fn encode_call(&self, site: &CalledSite) -> Result<vcf::variant::RecordBuf> {
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
        let allele_counts = site
            .allele_counts()
            .iter()
            .skip(1)
            .map(|&value| checked_info_integer(value, AC).map(Some))
            .collect::<Result<Vec<_>>>()?;
        let mut info = RecordInfo::default();
        if !allele_counts.is_empty() {
            info.insert(
                AC.to_owned(),
                Some(InfoValue::Array(InfoArray::Integer(allele_counts))),
            );
        }
        info.insert(
            AN.to_owned(),
            Some(InfoValue::Integer(
                i32::try_from(site.allele_number())
                    .map_err(|_| invalid("INFO/AN exceeds the VCF integer range"))?,
            )),
        );
        insert_indel_info(&self.header, &mut info, site.indel_summary())?;
        let include_dp = self.header.formats().contains_key(DP);
        let include_ad = self.header.formats().contains_key(AD);
        let include_qs = self.header.formats().contains_key(QUALITY_SUM);
        let include_gp = self.header.formats().contains_key(GP);
        let include_gq = self.header.formats().contains_key(GQ);
        let mut keys = vec![GT.to_owned(), PL.to_owned()];
        if include_dp {
            keys.push(DP.to_owned());
        }
        if include_ad {
            keys.push(AD.to_owned());
        }
        if include_qs {
            keys.push(QUALITY_SUM.to_owned());
        }
        if include_gp {
            keys.push(GP.to_owned());
        }
        if include_gq {
            keys.push(GQ.to_owned());
        }
        let keys = Keys::from_iter(keys);
        let values = site
            .samples()
            .iter()
            .map(|sample| {
                let genotype = sample.genotype().map(|alleles| {
                    Genotype::from_iter(
                        alleles
                            .iter()
                            .copied()
                            .map(|allele| GenotypeAllele::new(Some(allele), Phasing::Unphased)),
                    )
                });
                let evidence = sample.evidence();
                let mut values = vec![
                    genotype.map(SampleValue::Genotype),
                    sample
                        .phred_likelihoods()
                        .map(|values| checked_array(values, PL))
                        .transpose()?,
                ];
                if include_dp {
                    values.push(Some(SampleValue::Integer(checked_integer(
                        evidence.depth(),
                        DP,
                    )?)));
                }
                if include_ad {
                    values.push(Some(checked_array(evidence.allele_depths(), AD)?));
                }
                if include_qs {
                    values.push(Some(checked_array(
                        evidence.allele_quality_sums(),
                        QUALITY_SUM,
                    )?));
                }
                if include_gp {
                    values.push(sample.genotype_probabilities().map(|values| {
                        SampleValue::Array(SampleArray::Float(
                            values.iter().copied().map(Some).collect(),
                        ))
                    }));
                }
                if include_gq {
                    values.push(
                        sample
                            .genotype_quality()
                            .map(|value| SampleValue::Integer(i32::from(value))),
                    );
                }
                Ok(values)
            })
            .collect::<Result<Vec<_>>>()?;
        let mut builder = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(reference_name)
            .set_variant_start(position)
            .set_reference_bases(allele_text(site.reference())?)
            .set_alternate_bases(alternate_bases)
            .set_info(info)
            .set_samples(Samples::new(keys, values));
        if let Some(quality) = site.quality() {
            builder = builder.set_quality_score(quality);
        }
        Ok(builder.build())
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

fn require_optional_info(
    header: &vcf::Header,
    key: &str,
    number: InfoNumber,
    ty: InfoType,
) -> Result<()> {
    if header.infos().contains_key(key) {
        require_info(header, key, number, ty)?;
    }
    Ok(())
}

fn insert_indel_info(
    header: &vcf::Header,
    info: &mut RecordInfo,
    summary: Option<IndelSummary>,
) -> Result<()> {
    let Some(summary) = summary else {
        return Ok(());
    };
    require_info(header, INDEL, InfoNumber::Count(0), InfoType::Flag)?;
    require_info(header, IDV, InfoNumber::Count(1), InfoType::Integer)?;
    require_info(header, IMF, InfoNumber::Count(1), InfoType::Float)?;
    info.insert(INDEL.to_owned(), Some(InfoValue::Flag));
    info.insert(
        IDV.to_owned(),
        Some(InfoValue::Integer(checked_info_integer(
            summary.maximum_support(),
            IDV,
        )?)),
    );
    info.insert(
        IMF.to_owned(),
        Some(InfoValue::Float(summary.maximum_fraction())),
    );
    Ok(())
}

fn decode_indel_summary(record: &vcf::variant::RecordBuf) -> Result<Option<IndelSummary>> {
    let is_indel = match record.info().get(INDEL) {
        None => false,
        Some(Some(InfoValue::Flag)) => true,
        Some(_) => return Err(invalid("INFO/INDEL is not a flag")),
    };
    if !is_indel {
        if record.info().get(IDV).is_some() || record.info().get(IMF).is_some() {
            return Err(invalid(
                "INFO/IDV or INFO/IMF is present without INFO/INDEL",
            ));
        }
        return Ok(None);
    }
    let support = match record.info().get(IDV) {
        Some(Some(InfoValue::Integer(value))) => {
            u32::try_from(*value).map_err(|_| invalid("INFO/IDV contains a negative integer"))?
        }
        _ => return Err(invalid("indel record has no integer INFO/IDV value")),
    };
    let fraction = match record.info().get(IMF) {
        Some(Some(InfoValue::Float(value))) => *value,
        _ => return Err(invalid("indel record has no float INFO/IMF value")),
    };
    IndelSummary::new(support, fraction).map(Some)
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

fn checked_info_integer(value: u32, key: &str) -> Result<i32> {
    i32::try_from(value).map_err(|_| invalid(format!("INFO/{key} exceeds the VCF integer range")))
}

fn allele_text(allele: &Allele) -> Result<&str> {
    std::str::from_utf8(allele.as_bytes()).map_err(|_| invalid("allele is not UTF-8"))
}

fn invalid(message: impl Into<String>) -> CallError {
    CallError::InvalidLikelihoodVariant(message.into())
}

#[cfg(test)]
mod tests;
