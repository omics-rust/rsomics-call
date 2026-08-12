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
    Allele, CallError, CalledAnnotations, CalledSite, IndelSummary, LikelihoodSite, Ploidy,
    PriorAlleleCounts, Result, SampleAnnotations, SampleEvidence, SampleLikelihood,
    SiteAnnotations, model::SiteAnnotationValues,
};

const QUALITY_SUM: &str = "QS";
const INDEL: &str = "INDEL";
const IDV: &str = "IDV";
const IMF: &str = "IMF";
const VDB: &str = "VDB";
const RPBZ: &str = "RPBZ";
const MQBZ: &str = "MQBZ";
const BQBZ: &str = "BQBZ";
const MQSBZ: &str = "MQSBZ";
const NM: &str = "NM";
const NMBZ: &str = "NMBZ";
const SCBZ: &str = "SCBZ";
const FS: &str = "FS";
const SGB: &str = "SGB";
const MQ0F: &str = "MQ0F";
const I16: &str = "I16";
const ADF: &str = "ADF";
const ADR: &str = "ADR";
const SCR: &str = "SCR";
const SP: &str = "SP";
const QM: &str = "QM";
const PL: &str = "PL";
const DP: &str = "DP";
const AD: &str = "AD";
const AC: &str = "AC";
const AN: &str = "AN";
const GT: &str = "GT";
const GQ: &str = "GQ";
const GP: &str = "GP";
const DP4: &str = "DP4";
const MQ: &str = "MQ";
const PV4: &str = "PV4";
const END: &str = "END";
const MIN_DP: &str = "MIN_DP";

#[derive(Clone, Debug)]
pub struct LikelihoodVcfSchema {
    header: vcf::Header,
    prior_frequency_tags: Option<PriorFrequencyTags>,
}

#[derive(Clone, Debug)]
pub struct CalledVcfSchema {
    header: vcf::Header,
    prior_frequency_tags: Option<PriorFrequencyTags>,
}

#[derive(Clone, Debug)]
struct PriorFrequencyTags {
    total: String,
    alternates: String,
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
            .add_info(
                DP,
                Map::<Info>::new(InfoNumber::Count(1), InfoType::Integer, "Raw read depth"),
            )
            .add_info(
                VDB,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",
                ),
            )
            .add_info(
                RPBZ,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)",
                ),
            )
            .add_info(
                MQBZ,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)",
                ),
            )
            .add_info(
                BQBZ,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)",
                ),
            )
            .add_info(
                MQSBZ,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)",
                ),
            )
            .add_info(
                NM,
                Map::<Info>::new(
                    InfoNumber::Count(2),
                    InfoType::Float,
                    "Average number of mismatches in ref and alt reads",
                ),
            )
            .add_info(
                NMBZ,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better)",
                ),
            )
            .add_info(
                SCBZ,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)",
                ),
            )
            .add_info(
                FS,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Fisher's exact test P-value to detect strand bias",
                ),
            )
            .add_info(
                SGB,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Segregation based metric",
                ),
            )
            .add_info(
                MQ0F,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Float,
                    "Fraction of MQ0 reads (smaller is better)",
                ),
            )
            .add_info(
                I16,
                Map::<Info>::new(
                    InfoNumber::Count(16),
                    InfoType::Float,
                    "Auxiliary tag used for calling",
                ),
            )
            .add_info(
                AD,
                Map::<Info>::new(
                    InfoNumber::ReferenceAlternateBases,
                    InfoType::Integer,
                    "Total allelic depths (high-quality bases)",
                ),
            )
            .add_info(
                ADF,
                Map::<Info>::new(
                    InfoNumber::ReferenceAlternateBases,
                    InfoType::Integer,
                    "Total allelic depths on the forward strand (high-quality bases)",
                ),
            )
            .add_info(
                ADR,
                Map::<Info>::new(
                    InfoNumber::ReferenceAlternateBases,
                    InfoType::Integer,
                    "Total allelic depths on the reverse strand (high-quality bases)",
                ),
            )
            .add_info(
                SCR,
                Map::<Info>::new(
                    InfoNumber::Count(1),
                    InfoType::Integer,
                    "Number of soft-clipped reads",
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
                SP,
                Map::<Format>::new(
                    FormatNumber::Count(1),
                    FormatType::Integer,
                    "Phred-scaled strand bias P-value",
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
                ADF,
                Map::<Format>::new(
                    FormatNumber::ReferenceAlternateBases,
                    FormatType::Integer,
                    "Allelic depths on the forward strand (high-quality bases)",
                ),
            )
            .add_format(
                ADR,
                Map::<Format>::new(
                    FormatNumber::ReferenceAlternateBases,
                    FormatType::Integer,
                    "Allelic depths on the reverse strand (high-quality bases)",
                ),
            )
            .add_format(
                SCR,
                Map::<Format>::new(
                    FormatNumber::Count(1),
                    FormatType::Integer,
                    "Per-sample number of soft-clipped reads",
                ),
            )
            .add_format(
                QUALITY_SUM,
                Map::<Format>::new(
                    FormatNumber::ReferenceAlternateBases,
                    FormatType::Integer,
                    "Phred-score allele quality sum used by calling",
                ),
            )
            .add_format(
                QM,
                Map::<Format>::new(
                    FormatNumber::ReferenceAlternateBases,
                    FormatType::Integer,
                    "Phred-score allele quality mean",
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
        for (key, number, ty) in [
            (DP, InfoNumber::Count(1), InfoType::Integer),
            (VDB, InfoNumber::Count(1), InfoType::Float),
            (RPBZ, InfoNumber::Count(1), InfoType::Float),
            (MQBZ, InfoNumber::Count(1), InfoType::Float),
            (BQBZ, InfoNumber::Count(1), InfoType::Float),
            (MQSBZ, InfoNumber::Count(1), InfoType::Float),
            (NM, InfoNumber::Count(2), InfoType::Float),
            (NMBZ, InfoNumber::Count(1), InfoType::Float),
            (SCBZ, InfoNumber::Count(1), InfoType::Float),
            (FS, InfoNumber::Count(1), InfoType::Float),
            (SGB, InfoNumber::Count(1), InfoType::Float),
            (MQ0F, InfoNumber::Count(1), InfoType::Float),
            (I16, InfoNumber::Count(16), InfoType::Float),
            (AD, InfoNumber::ReferenceAlternateBases, InfoType::Integer),
            (ADF, InfoNumber::ReferenceAlternateBases, InfoType::Integer),
            (ADR, InfoNumber::ReferenceAlternateBases, InfoType::Integer),
            (SCR, InfoNumber::Count(1), InfoType::Integer),
        ] {
            require_optional_info(&header, key, number, ty)?;
        }
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
        for (key, number) in [
            (SP, FormatNumber::Count(1)),
            (ADF, FormatNumber::ReferenceAlternateBases),
            (ADR, FormatNumber::ReferenceAlternateBases),
            (SCR, FormatNumber::Count(1)),
            (QM, FormatNumber::ReferenceAlternateBases),
        ] {
            require_format(&header, key, number, FormatType::Integer, false)?;
        }
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
        Ok(Self {
            header,
            prior_frequency_tags: None,
        })
    }

    pub fn header(&self) -> &vcf::Header {
        &self.header
    }

    pub(crate) fn set_prior_frequency_tags(
        &mut self,
        total: impl Into<String>,
        alternates: impl Into<String>,
    ) -> Result<()> {
        let total = total.into();
        let alternates = alternates.into();
        if total.is_empty() || alternates.is_empty() || total == alternates {
            return Err(invalid(
                "prior-frequency INFO tags must be distinct and nonempty",
            ));
        }
        require_info(
            &self.header,
            &total,
            InfoNumber::Count(1),
            InfoType::Integer,
        )?;
        require_info(
            &self.header,
            &alternates,
            InfoNumber::AlternateBases,
            InfoType::Integer,
        )?;
        self.prior_frequency_tags = Some(PriorFrequencyTags { total, alternates });
        Ok(())
    }

    pub(crate) fn prior_frequency_tags(&self) -> Option<(&str, &str)> {
        self.prior_frequency_tags
            .as_ref()
            .map(|tags| (tags.total.as_str(), tags.alternates.as_str()))
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
        insert_prior_allele_counts(
            &mut info,
            self.prior_frequency_tags.as_ref(),
            site.prior_allele_counts(),
            false,
        )?;
        insert_indel_info(&self.header, &mut info, site.indel_summary())?;
        if let Some(annotations) = site.annotations() {
            insert_site_annotations(&self.header, &mut info, annotations)?;
        }
        let complete_sample_annotations = site
            .samples()
            .iter()
            .all(|sample| sample.evidence().annotations().is_some());
        if complete_sample_annotations {
            insert_total_depth_annotations(
                &self.header,
                &mut info,
                site.samples().iter().map(|sample| sample.evidence()),
            )?;
        }
        let mut keys = vec![
            PL.to_owned(),
            DP.to_owned(),
            AD.to_owned(),
            QUALITY_SUM.to_owned(),
        ];
        if complete_sample_annotations {
            keys = vec![
                PL.to_owned(),
                DP.to_owned(),
                SP.to_owned(),
                ADF.to_owned(),
                ADR.to_owned(),
                AD.to_owned(),
                SCR.to_owned(),
                QUALITY_SUM.to_owned(),
                QM.to_owned(),
            ];
        }
        let keys = Keys::from_iter(keys);
        let values = site
            .samples()
            .iter()
            .map(|sample| {
                let evidence = sample.evidence();
                let mut values = vec![
                    sample
                        .phred_likelihoods()
                        .map(|values| checked_array(values, PL))
                        .transpose()?,
                    Some(SampleValue::Integer(checked_integer(evidence.depth(), DP)?)),
                ];
                if let Some(annotations) = evidence.annotations() {
                    values.extend([
                        Some(SampleValue::Integer(checked_integer(
                            annotations.strand_bias(),
                            SP,
                        )?)),
                        Some(checked_array(annotations.forward_allele_depths(), ADF)?),
                        Some(checked_array(annotations.reverse_allele_depths(), ADR)?),
                        Some(checked_array(evidence.allele_depths(), AD)?),
                        Some(SampleValue::Integer(checked_integer(
                            annotations.soft_clipped_reads(),
                            SCR,
                        )?)),
                        Some(checked_array(evidence.allele_quality_sums(), QUALITY_SUM)?),
                        Some(checked_array(annotations.allele_quality_means(), QM)?),
                    ]);
                } else {
                    values.extend([
                        Some(checked_array(evidence.allele_depths(), AD)?),
                        Some(checked_array(evidence.allele_quality_sums(), QUALITY_SUM)?),
                    ]);
                }
                Ok(values)
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
        self.decode_likelihood_samples(record, 0..self.header.sample_names().len())
    }

    pub(crate) fn decode_selected_likelihood(
        &self,
        record: &vcf::variant::RecordBuf,
        sample_indices: &[usize],
    ) -> Result<LikelihoodSite> {
        self.decode_likelihood_samples(record, sample_indices.iter().copied())
    }

    fn decode_likelihood_samples(
        &self,
        record: &vcf::variant::RecordBuf,
        sample_indices: impl IntoIterator<Item = usize>,
    ) -> Result<LikelihoodSite> {
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

        let samples = sample_indices
            .into_iter()
            .map(|index| decode_sample(record.samples(), index, allele_count))
            .collect::<Result<Vec<_>>>()?;
        let mut site = LikelihoodSite::new(
            reference_sequence_id,
            position,
            reference,
            alternates,
            allele_quality_sums,
            samples,
        )?;
        if let Some(summary) = decode_indel_summary(record)? {
            site = site.with_indel_summary(summary);
        }
        if let Some(annotations) = decode_site_annotations(record)? {
            site = site.with_annotations(annotations);
        }
        if let Some(tags) = &self.prior_frequency_tags
            && let Some(counts) = decode_prior_allele_counts(record, tags, allele_count)?
        {
            site = site.with_prior_allele_counts(counts)?;
        }
        Ok(site)
    }
}

impl CalledVcfSchema {
    pub fn from_likelihood(schema: &LikelihoodVcfSchema) -> Self {
        Self::from_likelihood_inner(schema, true)
    }

    pub fn from_consensus_likelihood(schema: &LikelihoodVcfSchema) -> Self {
        Self::from_likelihood_inner(schema, false)
    }

    pub fn with_gvcf(mut self) -> Self {
        self.header.infos_mut().insert(
            END.to_owned(),
            Map::<Info>::new(
                InfoNumber::Count(1),
                InfoType::Integer,
                "End position of the variant described in this record",
            ),
        );
        self.header.infos_mut().insert(
            MIN_DP.to_owned(),
            Map::<Info>::new(
                InfoNumber::Count(1),
                InfoType::Integer,
                "Minimum per-sample depth in this gVCF block",
            ),
        );
        self
    }

    fn from_likelihood_inner(schema: &LikelihoodVcfSchema, include_probabilities: bool) -> Self {
        let mut header = schema.header.clone();
        header.infos_mut().shift_remove(QUALITY_SUM);
        header.infos_mut().shift_remove(I16);
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
        header.infos_mut().insert(
            DP4.to_owned(),
            Map::<Info>::new(
                InfoNumber::Count(4),
                InfoType::Integer,
                "High-quality ref-forward, ref-reverse, alt-forward and alt-reverse bases",
            ),
        );
        header.infos_mut().insert(
            MQ.to_owned(),
            Map::<Info>::new(
                InfoNumber::Count(1),
                InfoType::Integer,
                if include_probabilities {
                    "Average mapping quality"
                } else {
                    "Root-mean-square mapping quality of covering reads"
                },
            ),
        );
        header.infos_mut().insert(
            PV4.to_owned(),
            Map::<Info>::new(
                InfoNumber::Count(4),
                InfoType::Float,
                "P-values for strand, base-quality, mapping-quality and tail-distance bias",
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
        Self {
            header,
            prior_frequency_tags: schema.prior_frequency_tags.clone(),
        }
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
        if let Some(gvcf) = site.gvcf()
            && gvcf.is_collapsed()
        {
            return encode_gvcf_block(
                &self.header,
                reference_name,
                position,
                alternate_bases,
                site,
                gvcf,
            );
        }
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
        insert_prior_allele_counts(
            &mut info,
            self.prior_frequency_tags.as_ref(),
            site.prior_allele_counts(),
            true,
        )?;
        insert_indel_info(&self.header, &mut info, site.indel_summary())?;
        if let Some(annotations) = site.annotations() {
            insert_called_annotations(&self.header, &mut info, annotations)?;
        }
        if let Some(gvcf) = site.gvcf() {
            insert_gvcf_info(&self.header, &mut info, gvcf)?;
        }
        let complete_sample_annotations = site
            .samples()
            .iter()
            .all(|sample| sample.evidence().annotations().is_some());
        if complete_sample_annotations {
            insert_total_depth_annotations(
                &self.header,
                &mut info,
                site.samples().iter().map(|sample| sample.evidence()),
            )?;
        }
        let include_dp = self.header.formats().contains_key(DP);
        let include_ad = self.header.formats().contains_key(AD);
        let include_qs = self.header.formats().contains_key(QUALITY_SUM);
        let include_sample_annotations = complete_sample_annotations
            && [SP, ADF, ADR, SCR, QM]
                .iter()
                .all(|key| self.header.formats().contains_key(*key));
        let include_gp = self.header.formats().contains_key(GP);
        let include_gq = self.header.formats().contains_key(GQ);
        let mut keys = vec![GT.to_owned(), PL.to_owned()];
        if include_dp {
            keys.push(DP.to_owned());
        }
        if include_sample_annotations {
            keys.extend([SP, ADF, ADR].map(str::to_owned));
        }
        if include_ad {
            keys.push(AD.to_owned());
        }
        if include_sample_annotations {
            keys.push(SCR.to_owned());
        }
        if include_qs {
            keys.push(QUALITY_SUM.to_owned());
        }
        if include_sample_annotations {
            keys.push(QM.to_owned());
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
                if include_sample_annotations {
                    let annotations = evidence.annotations().unwrap();
                    values.extend([
                        Some(SampleValue::Integer(checked_integer(
                            annotations.strand_bias(),
                            SP,
                        )?)),
                        Some(checked_array(annotations.forward_allele_depths(), ADF)?),
                        Some(checked_array(annotations.reverse_allele_depths(), ADR)?),
                    ]);
                }
                if include_ad {
                    values.push(Some(checked_array(evidence.allele_depths(), AD)?));
                }
                if include_sample_annotations {
                    values.push(Some(SampleValue::Integer(checked_integer(
                        evidence.annotations().unwrap().soft_clipped_reads(),
                        SCR,
                    )?)));
                }
                if include_qs {
                    values.push(Some(checked_array(
                        evidence.allele_quality_sums(),
                        QUALITY_SUM,
                    )?));
                }
                if include_sample_annotations {
                    values.push(Some(checked_array(
                        evidence.annotations().unwrap().allele_quality_means(),
                        QM,
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

fn encode_gvcf_block(
    header: &vcf::Header,
    reference_name: &str,
    position: Position,
    alternate_bases: AlternateBases,
    site: &CalledSite,
    gvcf: crate::GvcfSite,
) -> Result<vcf::variant::RecordBuf> {
    require_format(header, GT, FormatNumber::Count(1), FormatType::String, true)?;
    require_format(
        header,
        DP,
        FormatNumber::Count(1),
        FormatType::Integer,
        true,
    )?;
    let include_pl = site
        .samples()
        .iter()
        .any(|sample| sample.phred_likelihoods().is_some());
    if include_pl {
        require_format(header, PL, FormatNumber::Samples, FormatType::Integer, true)?;
    }

    let mut info = RecordInfo::default();
    insert_gvcf_info(header, &mut info, gvcf)?;
    let mut keys = vec![GT.to_owned()];
    if include_pl {
        keys.push(PL.to_owned());
    }
    keys.push(DP.to_owned());
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
            let mut values = vec![genotype.map(SampleValue::Genotype)];
            if include_pl {
                values.push(
                    sample
                        .phred_likelihoods()
                        .map(|values| checked_array(values, PL))
                        .transpose()?,
                );
            }
            values.push(Some(SampleValue::Integer(checked_integer(
                sample.evidence().depth(),
                DP,
            )?)));
            Ok(values)
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(vcf::variant::RecordBuf::builder()
        .set_reference_sequence_name(reference_name)
        .set_variant_start(position)
        .set_reference_bases(allele_text(site.reference())?)
        .set_alternate_bases(alternate_bases)
        .set_info(info)
        .set_samples(Samples::new(Keys::from_iter(keys), values))
        .build())
}

fn insert_gvcf_info(
    header: &vcf::Header,
    info: &mut RecordInfo,
    gvcf: crate::GvcfSite,
) -> Result<()> {
    if let Some(end_position) = gvcf.end_position() {
        require_info(header, END, InfoNumber::Count(1), InfoType::Integer)?;
        let end = end_position
            .checked_add(1)
            .and_then(|value| i32::try_from(value).ok())
            .ok_or_else(|| invalid("gVCF block end exceeds the VCF integer range"))?;
        info.insert(END.to_owned(), Some(InfoValue::Integer(end)));
    }
    require_info(header, MIN_DP, InfoNumber::Count(1), InfoType::Integer)?;
    info.insert(
        MIN_DP.to_owned(),
        Some(InfoValue::Integer(checked_info_integer(
            gvcf.minimum_depth(),
            MIN_DP,
        )?)),
    );
    Ok(())
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

fn decode_prior_allele_counts(
    record: &vcf::variant::RecordBuf,
    tags: &PriorFrequencyTags,
    allele_count: usize,
) -> Result<Option<PriorAlleleCounts>> {
    let Some(total) = info_integer(record, &tags.total)? else {
        return Ok(None);
    };
    let Some(value) = record.info().get(&tags.alternates) else {
        return Ok(None);
    };
    let Some(InfoValue::Array(InfoArray::Integer(values))) = value else {
        return Err(invalid(format!(
            "INFO/{} is not an integer array",
            tags.alternates
        )));
    };
    if values.len() != allele_count - 1 {
        return Err(invalid(format!(
            "INFO/{} does not match the record alternate alleles",
            tags.alternates
        )));
    }
    let alternates = values
        .iter()
        .map(|value| {
            let value = value.ok_or_else(|| {
                invalid(format!("INFO/{} contains a missing value", tags.alternates))
            })?;
            u32::try_from(value).map_err(|_| {
                invalid(format!(
                    "INFO/{} contains a negative integer",
                    tags.alternates
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;
    PriorAlleleCounts::new(total, alternates).map(Some)
}

fn insert_prior_allele_counts(
    info: &mut RecordInfo,
    tags: Option<&PriorFrequencyTags>,
    counts: Option<&PriorAlleleCounts>,
    called: bool,
) -> Result<()> {
    let (Some(tags), Some(counts)) = (tags, counts) else {
        return Ok(());
    };
    if !called || tags.total != AN {
        info.insert(
            tags.total.clone(),
            Some(InfoValue::Integer(checked_info_integer(
                counts.total(),
                &tags.total,
            )?)),
        );
    }
    if (!called || tags.alternates != AC) && !counts.alternates().is_empty() {
        let values = counts
            .alternates()
            .iter()
            .map(|&count| checked_info_integer(count, &tags.alternates).map(Some))
            .collect::<Result<Vec<_>>>()?;
        info.insert(
            tags.alternates.clone(),
            Some(InfoValue::Array(InfoArray::Integer(values))),
        );
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

fn insert_site_annotations(
    header: &vcf::Header,
    info: &mut RecordInfo,
    annotations: &SiteAnnotations,
) -> Result<()> {
    insert_info_integer(header, info, DP, annotations.raw_depth())?;
    insert_info_float_array(header, info, I16, annotations.auxiliary())?;
    insert_optional_info_float(header, info, VDB, annotations.variant_distance_bias())?;
    insert_optional_info_float(header, info, SGB, annotations.segregation_bias())?;
    if let Some(values) = annotations.average_mismatches() {
        insert_info_float_array(header, info, NM, &values)?;
    }
    insert_optional_info_float(header, info, RPBZ, annotations.read_position_bias())?;
    insert_optional_info_float(header, info, MQBZ, annotations.mapping_quality_bias())?;
    insert_optional_info_float(
        header,
        info,
        MQSBZ,
        annotations.mapping_quality_strand_bias(),
    )?;
    insert_optional_info_float(header, info, BQBZ, annotations.base_quality_bias())?;
    insert_optional_info_float(header, info, NMBZ, annotations.mismatch_bias())?;
    insert_optional_info_float(header, info, SCBZ, annotations.soft_clip_bias())?;
    insert_optional_info_float(header, info, FS, annotations.strand_bias())?;
    insert_info_float(
        header,
        info,
        MQ0F,
        annotations.zero_mapping_quality_fraction(),
    )
}

fn insert_called_annotations(
    header: &vcf::Header,
    info: &mut RecordInfo,
    annotations: &CalledAnnotations,
) -> Result<()> {
    insert_site_annotations(header, info, annotations.pileup())?;
    let strand_depths = annotations
        .strand_depths()
        .iter()
        .copied()
        .map(u64::from)
        .collect::<Vec<_>>();
    insert_info_integer_array(header, info, DP4, &strand_depths)?;
    insert_info_integer(header, info, MQ, annotations.mapping_quality())?;
    if let Some(values) = annotations.bias_probabilities() {
        insert_info_float_array(header, info, PV4, values)?;
    }
    Ok(())
}

fn insert_total_depth_annotations<'a>(
    header: &vcf::Header,
    info: &mut RecordInfo,
    samples: impl IntoIterator<Item = &'a SampleEvidence>,
) -> Result<()> {
    let samples = samples.into_iter().collect::<Vec<_>>();
    let allele_count = samples[0].allele_depths().len();
    let mut forward = vec![0u64; allele_count];
    let mut reverse = vec![0u64; allele_count];
    let mut soft_clipped_reads = 0u64;
    for sample in samples {
        let annotations = sample.annotations().unwrap();
        for (total, &value) in forward.iter_mut().zip(annotations.forward_allele_depths()) {
            *total += u64::from(value);
        }
        for (total, &value) in reverse.iter_mut().zip(annotations.reverse_allele_depths()) {
            *total += u64::from(value);
        }
        soft_clipped_reads += u64::from(annotations.soft_clipped_reads());
    }
    let depths = forward
        .iter()
        .zip(&reverse)
        .map(|(&forward, &reverse)| forward + reverse)
        .collect::<Vec<_>>();
    insert_info_integer_array(header, info, ADF, &forward)?;
    insert_info_integer_array(header, info, ADR, &reverse)?;
    insert_info_integer_array(header, info, AD, &depths)?;
    insert_info_integer(
        header,
        info,
        SCR,
        u32::try_from(soft_clipped_reads)
            .map_err(|_| invalid("INFO/SCR exceeds the VCF integer range"))?,
    )
}

fn insert_info_integer(
    header: &vcf::Header,
    info: &mut RecordInfo,
    key: &str,
    value: u32,
) -> Result<()> {
    if header.infos().contains_key(key) {
        info.insert(
            key.to_owned(),
            Some(InfoValue::Integer(checked_info_integer(value, key)?)),
        );
    }
    Ok(())
}

fn insert_info_integer_array(
    header: &vcf::Header,
    info: &mut RecordInfo,
    key: &str,
    values: &[u64],
) -> Result<()> {
    if header.infos().contains_key(key) {
        let values = values
            .iter()
            .map(|&value| {
                i32::try_from(value)
                    .map(Some)
                    .map_err(|_| invalid(format!("INFO/{key} exceeds the VCF integer range")))
            })
            .collect::<Result<Vec<_>>>()?;
        info.insert(
            key.to_owned(),
            Some(InfoValue::Array(InfoArray::Integer(values))),
        );
    }
    Ok(())
}

fn insert_info_float(
    header: &vcf::Header,
    info: &mut RecordInfo,
    key: &str,
    value: f32,
) -> Result<()> {
    if header.infos().contains_key(key) {
        info.insert(key.to_owned(), Some(InfoValue::Float(value)));
    }
    Ok(())
}

fn insert_optional_info_float(
    header: &vcf::Header,
    info: &mut RecordInfo,
    key: &str,
    value: Option<f32>,
) -> Result<()> {
    if let Some(value) = value {
        insert_info_float(header, info, key, value)?;
    }
    Ok(())
}

fn insert_info_float_array(
    header: &vcf::Header,
    info: &mut RecordInfo,
    key: &str,
    values: &[f32],
) -> Result<()> {
    if header.infos().contains_key(key) {
        info.insert(
            key.to_owned(),
            Some(InfoValue::Array(InfoArray::Float(
                values.iter().copied().map(Some).collect(),
            ))),
        );
    }
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
    let mut evidence = SampleEvidence::new(depth, allele_depths, allele_quality_sums)?;
    let forward = sample_integer_array(samples, ADF, index)?;
    let reverse = sample_integer_array(samples, ADR, index)?;
    match (forward, reverse) {
        (None, None) => {}
        (Some(forward), Some(reverse))
            if forward.len() == allele_count && reverse.len() == allele_count =>
        {
            let quality_means =
                sample_quality_mean_array(samples, index)?.unwrap_or_else(|| vec![0; allele_count]);
            if quality_means.len() != allele_count {
                return Err(invalid("FORMAT/QM does not match the record alleles"));
            }
            let annotations = SampleAnnotations::new(
                forward,
                reverse,
                quality_means,
                sample_integer(samples, SP, index)?.unwrap_or(0),
                sample_integer(samples, SCR, index)?.unwrap_or(0),
            )?;
            evidence = evidence.with_annotations(annotations)?;
        }
        _ => {
            return Err(invalid(
                "FORMAT/ADF and FORMAT/ADR must be present together and match the record alleles",
            ));
        }
    }
    SampleLikelihood::new(
        ploidy,
        phred_likelihoods.map(Vec::into_boxed_slice),
        evidence,
    )
}

fn decode_site_annotations(record: &vcf::variant::RecordBuf) -> Result<Option<SiteAnnotations>> {
    let Some(auxiliary) = info_float_array(record, I16)? else {
        return Ok(None);
    };
    let auxiliary: [f32; 16] = auxiliary
        .try_into()
        .map_err(|_| invalid("INFO/I16 must contain 16 values"))?;
    let average_mismatches = info_float_array(record, NM)?
        .map(|values| {
            values
                .try_into()
                .map_err(|_| invalid("INFO/NM must contain two values"))
        })
        .transpose()?;
    let raw_depth = info_integer(record, DP)?.unwrap_or_else(|| {
        auxiliary[..4]
            .iter()
            .copied()
            .map(|value| value.max(0.0) as u32)
            .sum()
    });
    SiteAnnotations::new(SiteAnnotationValues {
        raw_depth,
        auxiliary,
        variant_distance_bias: info_float(record, VDB)?,
        read_position_bias: info_float(record, RPBZ)?,
        mapping_quality_bias: info_float(record, MQBZ)?,
        base_quality_bias: info_float(record, BQBZ)?,
        mapping_quality_strand_bias: info_float(record, MQSBZ)?,
        mismatch_bias: info_float(record, NMBZ)?,
        soft_clip_bias: info_float(record, SCBZ)?,
        strand_bias: info_float(record, FS)?,
        segregation_bias: info_float(record, SGB)?,
        zero_mapping_quality_fraction: info_float(record, MQ0F)?.unwrap_or(0.0),
        average_mismatches,
    })
    .map(Some)
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

fn info_integer(record: &vcf::variant::RecordBuf, key: &str) -> Result<Option<u32>> {
    match record.info().get(key) {
        None | Some(None) => Ok(None),
        Some(Some(InfoValue::Integer(value))) => u32::try_from(*value)
            .map(Some)
            .map_err(|_| invalid(format!("INFO/{key} contains a negative integer"))),
        Some(Some(_)) => Err(invalid(format!("INFO/{key} is not an integer"))),
    }
}

fn info_float(record: &vcf::variant::RecordBuf, key: &str) -> Result<Option<f32>> {
    match record.info().get(key) {
        None | Some(None) => Ok(None),
        Some(Some(InfoValue::Float(value))) if value.is_finite() => Ok(Some(*value)),
        Some(Some(InfoValue::Float(_))) => Err(invalid(format!("INFO/{key} is not finite"))),
        Some(Some(_)) => Err(invalid(format!("INFO/{key} is not a float"))),
    }
}

fn info_float_array(record: &vcf::variant::RecordBuf, key: &str) -> Result<Option<Vec<f32>>> {
    match record.info().get(key) {
        None | Some(None) => Ok(None),
        Some(Some(InfoValue::Array(InfoArray::Float(values)))) => values
            .iter()
            .map(|value| match value {
                Some(value) if value.is_finite() => Ok(*value),
                Some(_) => Err(invalid(format!("INFO/{key} contains a non-finite value"))),
                None => Err(invalid(format!("INFO/{key} contains a missing value"))),
            })
            .collect::<Result<Vec<_>>>()
            .map(Some),
        Some(Some(_)) => Err(invalid(format!("INFO/{key} is not a float array"))),
    }
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

fn sample_quality_mean_array(samples: &Samples, index: usize) -> Result<Option<Vec<u32>>> {
    let Some(series) = samples.select(QM) else {
        return Ok(None);
    };
    match series.get(index) {
        Some(None) => Ok(None),
        Some(Some(SampleValue::Array(SampleArray::Integer(values)))) => values
            .iter()
            .map(|value| match value {
                None | Some(i32::MAX) => Ok(0),
                Some(value) => u32::try_from(*value)
                    .map_err(|_| invalid("FORMAT/QM contains a negative integer")),
            })
            .collect::<Result<Vec<_>>>()
            .map(Some),
        Some(Some(_)) => Err(invalid("FORMAT/QM is not an integer array")),
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
