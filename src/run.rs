use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::{
    core::{Position, Region},
    fasta,
};
use rsomics_pileup::{PileupEngine, PileupOptions};

use crate::{
    AlignmentInput, AlignmentSet, CallError, CalledSite, LikelihoodSite, Nucleotide,
    ReferenceSequence, Result, SampleMap, SampleSelection, SnpLikelihoodConfig, SnpSiteBuilder,
};

pub struct SnpLikelihoodRun {
    alignments: AlignmentSet,
    reference: ReferenceCache,
    pileup: PileupEngine,
    sites: SnpSiteBuilder,
}

impl SnpLikelihoodRun {
    pub fn open(
        inputs: impl IntoIterator<Item = AlignmentInput>,
        reference: impl AsRef<Path>,
        selection: SampleSelection,
        pileup_options: PileupOptions,
        likelihood_config: SnpLikelihoodConfig,
    ) -> Result<Self> {
        let reference = reference.as_ref();
        let alignments = AlignmentSet::open(inputs, Some(reference), selection)?;
        let lengths = alignments
            .reference_sequences()
            .iter()
            .map(ReferenceSequence::length);
        let pileup = PileupEngine::new(lengths, pileup_options);
        let sites = SnpSiteBuilder::new(alignments.samples().samples().len(), likelihood_config)?;
        let reference = ReferenceCache::open(reference, alignments.reference_sequences())?;
        Ok(Self {
            alignments,
            reference,
            pileup,
            sites,
        })
    }

    pub fn reference_sequences(&self) -> &[ReferenceSequence] {
        self.alignments.reference_sequences()
    }

    pub fn samples(&self) -> &SampleMap {
        self.alignments.samples()
    }

    pub fn run(mut self, mut emit: impl FnMut(LikelihoodSite) -> Result<()>) -> Result<()> {
        while let Some((source_id, record)) = self.alignments.next_record()? {
            self.pileup.push_with_source(source_id, record)?;
            drain_sites(
                &mut self.pileup,
                &mut self.reference,
                &mut self.sites,
                self.alignments.samples(),
                &mut emit,
            )?;
        }
        self.pileup.finish()?;
        drain_sites(
            &mut self.pileup,
            &mut self.reference,
            &mut self.sites,
            self.alignments.samples(),
            &mut emit,
        )
    }

    pub fn run_called(
        self,
        mut call: impl FnMut(&LikelihoodSite) -> Result<CalledSite>,
        mut emit: impl FnMut(CalledSite) -> Result<()>,
    ) -> Result<()> {
        self.run(|site| emit(call(&site)?))
    }
}

fn drain_sites(
    pileup: &mut PileupEngine,
    reference: &mut ReferenceCache,
    sites: &mut SnpSiteBuilder,
    samples: &SampleMap,
    emit: &mut impl FnMut(LikelihoodSite) -> Result<()>,
) -> Result<()> {
    pileup.drain(|column| {
        let reference_base = reference.base(column.reference_id(), column.position())?;
        let site = sites.build(column, reference_base, |source_id, record| {
            samples.sample_index(source_id, record)
        })?;
        emit(site)
    })
}

struct ReferenceCache {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    path: PathBuf,
    references: Box<[(Box<[u8]>, u64)]>,
    reference_id: Option<usize>,
    sequence_start: usize,
    sequence: Vec<u8>,
}

impl ReferenceCache {
    fn open(path: &Path, references: &[ReferenceSequence]) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|error| reference_error(path, error))?;
        Ok(Self {
            reader,
            path: path.to_path_buf(),
            references: references
                .iter()
                .map(|reference| {
                    (
                        reference.name().to_vec().into_boxed_slice(),
                        reference.length(),
                    )
                })
                .collect(),
            reference_id: None,
            sequence_start: 0,
            sequence: Vec::new(),
        })
    }

    fn base(&mut self, reference_id: i32, position: i64) -> Result<Nucleotide> {
        const CHUNK_SIZE: usize = 1024 * 1024;

        let reference_id =
            usize::try_from(reference_id).map_err(|error| reference_error(&self.path, error))?;
        let position =
            usize::try_from(position).map_err(|error| reference_error(&self.path, error))?;
        let (name, length) = self
            .references
            .get(reference_id)
            .map(|(name, length)| (name.to_vec(), *length))
            .ok_or_else(|| reference_error(&self.path, "reference ID is absent"))?;
        let length = usize::try_from(length).map_err(|error| reference_error(&self.path, error))?;
        if self.reference_id != Some(reference_id)
            || position < self.sequence_start
            || position >= self.sequence_start + self.sequence.len()
        {
            let start = position / CHUNK_SIZE * CHUNK_SIZE;
            let end = start
                .checked_add(CHUNK_SIZE)
                .map_or(length, |end| end.min(length));
            let interval_start = Position::try_from(start + 1)
                .map_err(|error| reference_error(&self.path, error))?;
            let interval_end =
                Position::try_from(end).map_err(|error| reference_error(&self.path, error))?;
            let record = self
                .reader
                .query(&Region::new(name, interval_start..=interval_end))
                .map_err(|error| reference_error(&self.path, error))?;
            self.sequence.clear();
            self.sequence.extend_from_slice(record.sequence().as_ref());
            self.reference_id = Some(reference_id);
            self.sequence_start = start;
        }
        let offset = position
            .checked_sub(self.sequence_start)
            .ok_or_else(|| reference_error(&self.path, "invalid reference cache position"))?;
        let base = *self.sequence.get(offset).ok_or_else(|| {
            reference_error(&self.path, "pileup position is outside the reference")
        })?;
        Ok(match base.to_ascii_uppercase() {
            b'A' => Nucleotide::A,
            b'C' => Nucleotide::C,
            b'G' => Nucleotide::G,
            b'T' => Nucleotide::T,
            _ => Nucleotide::N,
        })
    }
}

fn reference_error(path: &Path, error: impl std::fmt::Display) -> CallError {
    CallError::ReferenceInput {
        path: path.display().to_string(),
        message: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use tempfile::NamedTempFile;

    use super::*;

    fn sam_file(sample: &str, base: char) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:20\n\
             @RG\tID:rg\tSM:{sample}\n\
             read\t0\tchr1\t1\t60\t1M\t*\t0\t0\t{base}\tI\tRG:Z:rg"
        )
        .unwrap();
        file
    }

    fn repeated_sam_file(count: usize) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "@HD\tVN:1.6\tSO:coordinate\n\
             @SQ\tSN:chr1\tLN:20\n\
             @RG\tID:rg\tSM:S1"
        )
        .unwrap();
        for index in 0..count {
            writeln!(
                file,
                "read{index}\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\tRG:Z:rg"
            )
            .unwrap();
        }
        file
    }

    #[test]
    fn streams_cram_records_into_likelihood_sites() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let reference = fixtures.join("alignment-reference.fa");
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(
                3,
                fixtures.join("alignment.cram"),
                "input",
            )],
            reference,
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap();
        assert_eq!(
            run.samples()
                .samples()
                .iter()
                .map(AsRef::as_ref)
                .collect::<Vec<_>>(),
            ["S1"]
        );

        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 2);
        assert_eq!(
            (sites[0].position(), sites[0].reference().as_bytes()),
            (1, b"C".as_slice())
        );
        assert_eq!(
            (sites[1].position(), sites[1].reference().as_bytes()),
            (2, b"G".as_slice())
        );
        for site in sites {
            assert_eq!(site.samples()[0].phred_likelihoods(), Some(&[0, 3, 40][..]));
            assert_eq!(site.samples()[0].evidence().depth(), 1);
        }
    }

    #[test]
    fn matches_bcftools_1_24_multisample_snp_likelihoods() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let first = sam_file("S1", 'A');
        let second = sam_file("S2", 'G');
        let run = SnpLikelihoodRun::open(
            [
                AlignmentInput::new(1, first.path(), "first"),
                AlignmentInput::new(2, second.path(), "second"),
            ],
            fixtures.join("alignment-reference.fa"),
            SampleSelection::default(),
            PileupOptions::default(),
            SnpLikelihoodConfig::default(),
        )
        .unwrap();

        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 1);
        let site = &sites[0];
        assert_eq!(site.position(), 0);
        assert_eq!(site.reference().as_bytes(), b"A");
        assert_eq!(
            site.alternates()
                .iter()
                .map(|allele| allele.as_bytes())
                .collect::<Vec<_>>(),
            [b"G".as_slice(), b"<*>".as_slice()]
        );
        assert_eq!(
            site.samples()[0].phred_likelihoods(),
            Some(&[0, 3, 40, 3, 40, 40][..])
        );
        assert_eq!(
            site.samples()[1].phred_likelihoods(),
            Some(&[40, 3, 0, 40, 3, 40][..])
        );
    }

    #[test]
    fn fused_call_matches_materialized_typed_pipeline() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let reference = fixtures.join("alignment-reference.fa");
        let first = sam_file("S1", 'A');
        let second = sam_file("S2", 'G');
        let open = || {
            SnpLikelihoodRun::open(
                [
                    AlignmentInput::new(1, first.path(), "first"),
                    AlignmentInput::new(2, second.path(), "second"),
                ],
                &reference,
                SampleSelection::default(),
                PileupOptions::default(),
                SnpLikelihoodConfig::default(),
            )
            .unwrap()
        };
        let caller = crate::MultiallelicCaller::default();
        let mut materialized = Vec::new();
        open()
            .run(|site| {
                materialized.push(site);
                Ok(())
            })
            .unwrap();
        let expected = materialized
            .iter()
            .map(|site| caller.call(site).unwrap())
            .collect::<Vec<_>>();
        let mut fused = Vec::new();

        open()
            .run_called(
                |site| caller.call(site),
                |site| {
                    fused.push(site);
                    Ok(())
                },
            )
            .unwrap();

        assert_eq!(fused, expected);
    }

    #[test]
    fn matches_bcftools_1_24_per_input_depth_limit() {
        let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
        let input = repeated_sam_file(5);
        let run = SnpLikelihoodRun::open(
            [AlignmentInput::new(1, input.path(), "input")],
            fixtures.join("alignment-reference.fa"),
            SampleSelection::default(),
            PileupOptions {
                maximum_depth_per_source: Some(2),
                ..PileupOptions::default()
            },
            SnpLikelihoodConfig::default(),
        )
        .unwrap();

        let mut sites = Vec::new();
        run.run(|site| {
            sites.push(site);
            Ok(())
        })
        .unwrap();

        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].samples()[0].evidence().depth(), 2);
        assert_eq!(
            sites[0].samples()[0].phred_likelihoods(),
            Some(&[0, 6, 73][..])
        );
    }
}
