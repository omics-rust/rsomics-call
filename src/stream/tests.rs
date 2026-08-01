use std::io::{self, Write};

use super::*;
use crate::{Allele, MultiallelicCaller, Ploidy, SampleEvidence, SampleLikelihood};

const FORMATS: [VariantOutputFormat; 4] = [
    VariantOutputFormat::Vcf,
    VariantOutputFormat::VcfBgzf,
    VariantOutputFormat::BcfRaw,
    VariantOutputFormat::BcfBgzf,
];

#[test]
fn roundtrips_likelihood_streams_in_all_formats() {
    let (schema, site) = fixture();

    for format in FORMATS {
        let mut writer = LikelihoodVariantWriter::new(Vec::new(), schema.clone(), format).unwrap();
        writer.write_site(&site).unwrap();
        let data = writer.finish().unwrap();
        let mut reader = LikelihoodVariantReader::new(&data[..]).unwrap();

        assert_eq!(
            reader.schema().header().file_format(),
            schema.header().file_format()
        );
        assert_eq!(
            reader
                .schema()
                .header()
                .contigs()
                .keys()
                .collect::<Vec<_>>(),
            schema.header().contigs().keys().collect::<Vec<_>>()
        );
        assert_eq!(
            reader
                .schema()
                .header()
                .sample_names()
                .iter()
                .collect::<Vec<_>>(),
            schema.header().sample_names().iter().collect::<Vec<_>>()
        );
        assert_eq!(reader.read_site().unwrap(), Some(site.clone()));
        assert_eq!(reader.read_site().unwrap(), None);
    }
}

#[test]
fn writes_called_streams_in_all_formats() {
    let (likelihood_schema, site) = fixture();
    let called = MultiallelicCaller::default().call(&site).unwrap();

    for format in FORMATS {
        let schema = CalledVcfSchema::from_likelihood(&likelihood_schema);
        let mut writer = CalledVariantWriter::new(Vec::new(), schema, format).unwrap();
        writer.write_site(&called).unwrap();
        let data = writer.finish().unwrap();
        let mut reader = variant::io::Reader::new(&data[..]).unwrap();
        let header = reader.read_header().unwrap();
        let mut record = variant::Record::default();

        assert_ne!(reader.read_record(&mut record).unwrap(), 0);
        let record = vcf::variant::RecordBuf::try_from_variant_record(&header, &record).unwrap();
        assert_eq!(record.reference_sequence_name(), "chr1");
        assert!(record.info().get("AN").is_some());
        assert_eq!(
            reader.read_record(&mut variant::Record::default()).unwrap(),
            0
        );
    }
}

#[test]
fn rejects_truncated_compressed_input() {
    let (schema, site) = fixture();
    let mut writer =
        LikelihoodVariantWriter::new(Vec::new(), schema, VariantOutputFormat::BcfBgzf).unwrap();
    writer.write_site(&site).unwrap();
    let mut data = writer.finish().unwrap();
    data.truncate(data.len() / 2);

    let result = LikelihoodVariantReader::new(&data[..]).and_then(|mut reader| reader.read_site());
    assert!(matches!(
        result,
        Err(CallError::LikelihoodVariantInput(_) | CallError::LikelihoodVariantRecord { .. })
    ));
}

#[test]
fn reports_the_invalid_record_number() {
    let src = b"##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100>\n\
##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n\
chr1\t1\t.\tA\tG\t.\t.\tQS=1,1\tPL\t40,3,0\n\
chr1\t2\t.\tA\tG\t.\t.\t.\tPL\t40,3,0\n";
    let mut reader = LikelihoodVariantReader::new(&src[..]).unwrap();

    reader.read_site().unwrap();
    assert!(matches!(
        reader.read_site(),
        Err(CallError::LikelihoodVariantRecord { record: 2, .. })
    ));
}

#[test]
fn reports_deferred_output_failures() {
    let (schema, site) = fixture();
    let mut writer =
        LikelihoodVariantWriter::new(RejectWrites, schema, VariantOutputFormat::VcfBgzf).unwrap();
    writer.write_site(&site).unwrap();

    assert!(matches!(writer.finish(), Err(CallError::VariantOutput(_))));
}

#[test]
fn projects_likelihood_samples_during_decode() {
    let (schema, site) = multi_sample_fixture();
    for format in FORMATS {
        let mut writer = LikelihoodVariantWriter::new(Vec::new(), schema.clone(), format).unwrap();
        writer.write_site(&site).unwrap();
        let data = writer.finish().unwrap();

        let mut reader = LikelihoodVariantReader::new(&data[..])
            .unwrap()
            .select_samples(["third", "first"])
            .unwrap();
        assert_eq!(
            reader
                .schema()
                .header()
                .sample_names()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["third", "first"]
        );
        let projected = reader.read_site().unwrap().unwrap();
        assert_eq!(projected.samples().len(), 2);
        assert_eq!(
            projected.samples()[0].phred_likelihoods(),
            site.samples()[2].phred_likelihoods()
        );
        assert_eq!(
            projected.samples()[1].phred_likelihoods(),
            site.samples()[0].phred_likelihoods()
        );

        let mut reader = LikelihoodVariantReader::new(&data[..])
            .unwrap()
            .exclude_samples(["second"])
            .unwrap();
        assert_eq!(
            reader
                .schema()
                .header()
                .sample_names()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            ["first", "third"]
        );
        assert_eq!(reader.read_site().unwrap().unwrap().samples().len(), 2);
    }
}

#[test]
fn rejects_invalid_or_late_likelihood_sample_selection() {
    let (schema, site) = multi_sample_fixture();
    let mut writer =
        LikelihoodVariantWriter::new(Vec::new(), schema, VariantOutputFormat::Vcf).unwrap();
    writer.write_site(&site).unwrap();
    let data = writer.finish().unwrap();

    assert!(matches!(
        LikelihoodVariantReader::new(&data[..])
            .unwrap()
            .select_samples(["first", "first"]),
        Err(CallError::DuplicateSampleSelection(_))
    ));
    assert!(matches!(
        LikelihoodVariantReader::new(&data[..])
            .unwrap()
            .exclude_samples(["absent"]),
        Err(CallError::MissingSelectedSample(_))
    ));
    assert!(matches!(
        LikelihoodVariantReader::new(&data[..])
            .unwrap()
            .exclude_samples(["first", "second", "third"]),
        Err(CallError::InvalidSampleCount)
    ));
    let mut reader = LikelihoodVariantReader::new(&data[..]).unwrap();
    reader.read_site().unwrap();
    assert!(matches!(
        reader.select_samples(["first"]),
        Err(CallError::LateLikelihoodSampleSelection)
    ));
}

fn fixture() -> (LikelihoodVcfSchema, LikelihoodSite) {
    let schema = LikelihoodVcfSchema::new([(b"chr1".as_slice(), 100)], ["sample"]).unwrap();
    let site = LikelihoodSite::new(
        0,
        9,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"G"[..]).unwrap()],
        [1.0, 1.0],
        [SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            [40, 3, 0],
            SampleEvidence::new(1, [0, 1], [0, 40]).unwrap(),
        )
        .unwrap()],
    )
    .unwrap();
    (schema, site)
}

fn multi_sample_fixture() -> (LikelihoodVcfSchema, LikelihoodSite) {
    let schema =
        LikelihoodVcfSchema::new([(b"chr1".as_slice(), 100)], ["first", "second", "third"])
            .unwrap();
    let samples = [[0, 3, 40], [40, 3, 0], [20, 0, 20]].map(|likelihoods| {
        SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            likelihoods,
            SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
        )
        .unwrap()
    });
    let site = LikelihoodSite::new(
        0,
        9,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"G"[..]).unwrap()],
        [1.0, 1.0],
        samples,
    )
    .unwrap();
    (schema, site)
}

struct RejectWrites;

impl Write for RejectWrites {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        Err(io::Error::other("rejected write"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
