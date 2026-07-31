use std::path::Path;

use noodles::vcf::variant::io::Write as _;

use super::*;

#[test]
fn decodes_bcftools_1_24_likelihood_record() {
    let data = include_bytes!("../../tests/golden/bcftools-1.24-likelihood.vcf");
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
    let called_schema = CalledVcfSchema::from_likelihood(&schema);
    assert!(!called_schema.header().infos().contains_key(QUALITY_SUM));
    let called_record = called_schema.encode_call(&called).unwrap();
    let mut called_vcf = Vec::new();
    let mut called_writer = vcf::io::Writer::new(&mut called_vcf);
    called_writer.write_header(called_schema.header()).unwrap();
    called_writer
        .write_variant_record(called_schema.header(), &called_record)
        .unwrap();
    let called_line = std::str::from_utf8(&called_vcf)
        .unwrap()
        .lines()
        .last()
        .unwrap();
    let called_fields = called_line.split('\t').collect::<Vec<_>>();
    assert_eq!(&called_fields[..5], ["chr1", "1", ".", "A", "G,C"]);
    assert_eq!(called_fields[7], "AC=2,1;AN=6");
    assert_eq!(called_fields[8], "GT:PL:DP:AD:QS:GP:GQ");
    assert!(called_fields[9].starts_with("0/1:0,3,40,3,40,40:"));
    let mut called_vcf_reader = vcf::io::Reader::new(&called_vcf[..]);
    let called_vcf_header = called_vcf_reader.read_header().unwrap();
    let mut called_vcf_record = vcf::variant::RecordBuf::default();
    called_vcf_reader
        .read_record_buf(&called_vcf_header, &mut called_vcf_record)
        .unwrap();
    assert_eq!(called_vcf_record, called_record);
    let mut called_bcf = Vec::new();
    let mut called_bcf_writer = noodles::bcf::io::Writer::from(&mut called_bcf);
    called_bcf_writer
        .write_header(called_schema.header())
        .unwrap();
    called_bcf_writer
        .write_variant_record(called_schema.header(), &called_record)
        .unwrap();
    let mut called_bcf_reader = noodles::bcf::io::Reader::from(&called_bcf[..]);
    let called_bcf_header = called_bcf_reader.read_header().unwrap();
    let mut called_bcf_record = vcf::variant::RecordBuf::default();
    called_bcf_reader
        .read_record_buf(&called_bcf_header, &mut called_bcf_record)
        .unwrap();
    assert_eq!(called_bcf_record, called_record);
    let mut minimal_header = schema.header().clone();
    minimal_header.formats_mut().shift_remove(DP);
    minimal_header.formats_mut().shift_remove(AD);
    minimal_header.formats_mut().shift_remove(QUALITY_SUM);
    let minimal_input = LikelihoodVcfSchema::from_header(minimal_header).unwrap();
    let minimal_output = CalledVcfSchema::from_likelihood(&minimal_input);
    let minimal_record = minimal_output.encode_call(&called).unwrap();
    assert_eq!(
        minimal_record
            .samples()
            .keys()
            .as_ref()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        [GT, PL, GP, GQ]
    );

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
        "chr1\t1\t.\tA\tG,C,<*>\t.\t.\tQS=1,1,1,0;DP=3;I16=1,0,2,0,40,1600,80,3200,60,3600,120,7200,0,0,0,0;VDB=0.02;SGB=-0.94647;RPBZ=0;MQBZ=0;BQBZ=0;SCBZ=0;MQ0F=0\tPL:DP:AD:QS\t0,3,40,3,40,40,3,40,40,40:1:1,0,0,0:40,0,0,0\t40,40,40,3,3,0,40,40,3,40:1:0,0,1,0:0,0,40,0\t40,3,0,40,3,40,40,3,40,40:1:0,1,0,0:0,40,0,0"
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

    assert_eq!(
        schema.header().file_format(),
        vcf::header::FileFormat::new(4, 2)
    );
    assert_eq!(schema.header().contigs().len(), 1);
    assert_eq!(schema.header().sample_names().len(), 2);
    assert!(schema.header().infos().contains_key(QUALITY_SUM));
    assert!(schema.header().infos().contains_key(INDEL));
    assert!(schema.header().infos().contains_key(IDV));
    assert!(schema.header().infos().contains_key(IMF));
    assert!(schema.header().formats().contains_key(PL));
    assert!(schema.header().formats().contains_key(DP));
    assert!(schema.header().formats().contains_key(AD));
}

#[test]
fn roundtrips_complete_pileup_annotations() {
    let schema = LikelihoodVcfSchema::new([(b"MX".as_slice(), 9)], ["sample"]).unwrap();
    let sample_annotations =
        SampleAnnotations::new([4, 1, 0], [1, 4, 0], [27, 20, 0], 7, 2).unwrap();
    let evidence = SampleEvidence::new(10, [5, 5, 0], [180, 100, 0])
        .unwrap()
        .with_annotations(sample_annotations)
        .unwrap();
    let annotations = SiteAnnotations::new(SiteAnnotationValues {
        raw_depth: 10,
        auxiliary: [
            4.0, 1.0, 1.0, 4.0, 200.0, 8000.0, 100.0, 2000.0, 260.0, 14_800.0, 170.0, 6100.0, 20.0,
            80.0, 20.0, 80.0,
        ],
        variant_distance_bias: Some(0.001_870_952_4),
        read_position_bias: Some(0.0),
        mapping_quality_bias: Some(-1.671_258_1),
        base_quality_bias: Some(-3.0),
        mapping_quality_strand_bias: Some(-2.785_43),
        mismatch_bias: Some(1.671_258_1),
        soft_clip_bias: Some(0.0),
        strand_bias: Some(0.099_206_35),
        segregation_bias: Some(-0.590_764_64),
        zero_mapping_quality_fraction: 0.0,
        average_mismatches: Some([0.6, 2.0]),
    })
    .unwrap();
    let site = LikelihoodSite::new(
        0,
        4,
        Allele::new(&b"A"[..]).unwrap(),
        [
            Allele::new(&b"C"[..]).unwrap(),
            Allele::new(&b"<*>"[..]).unwrap(),
        ],
        [0.642_857_13, 0.357_142_87, 0.0],
        [SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            [56, 0, 118, 71, 133, 180],
            evidence,
        )
        .unwrap()],
    )
    .unwrap()
    .with_annotations(annotations);

    let record = schema.encode_likelihood(&site).unwrap();
    assert_eq!(record.info().get(SCR), Some(Some(&InfoValue::Integer(2))));
    assert_eq!(
        record
            .samples()
            .keys()
            .as_ref()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        [PL, DP, SP, ADF, ADR, AD, SCR, QUALITY_SUM, QM]
    );
    assert_eq!(schema.decode_likelihood(&record).unwrap(), site);
}

#[test]
fn preserves_indel_summary_in_likelihood_and_called_records() {
    let schema = LikelihoodVcfSchema::new([(b"MX".as_slice(), 11)], ["sample"]).unwrap();
    let site = LikelihoodSite::new(
        0,
        4,
        Allele::new(&b"T"[..]).unwrap(),
        [Allele::new(&b"TC"[..]).unwrap()],
        [0.0, 1.0],
        [SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            [56, 6, 0],
            SampleEvidence::new(2, [0, 2], [0, 62]).unwrap(),
        )
        .unwrap()],
    )
    .unwrap()
    .with_indel_summary(IndelSummary::new(2, 1.0).unwrap());

    let record = schema.encode_likelihood(&site).unwrap();
    assert_eq!(record.info().get(INDEL), Some(Some(&InfoValue::Flag)));
    assert_eq!(record.info().get(IDV), Some(Some(&InfoValue::Integer(2))));
    assert_eq!(record.info().get(IMF), Some(Some(&InfoValue::Float(1.0))));
    assert_eq!(schema.decode_likelihood(&record).unwrap(), site);

    let mut data = Vec::new();
    let mut writer = noodles::bcf::io::Writer::from(&mut data);
    writer.write_header(schema.header()).unwrap();
    writer
        .write_variant_record(schema.header(), &record)
        .unwrap();
    let mut reader = noodles::bcf::io::Reader::from(&data[..]);
    let header = reader.read_header().unwrap();
    let decoded_schema = LikelihoodVcfSchema::from_header(header.clone()).unwrap();
    let mut decoded = vcf::variant::RecordBuf::default();
    reader.read_record_buf(&header, &mut decoded).unwrap();
    assert_eq!(decoded_schema.decode_likelihood(&decoded).unwrap(), site);

    let called = crate::MultiallelicCaller::default().call(&site).unwrap();
    let called_schema = CalledVcfSchema::from_likelihood(&schema);
    let called_record = called_schema.encode_call(&called).unwrap();
    assert_eq!(
        called_record.info().get(INDEL),
        Some(Some(&InfoValue::Flag))
    );
    assert_eq!(
        called_record.info().get(IDV),
        Some(Some(&InfoValue::Integer(2)))
    );
    assert_eq!(
        called_record.info().get(IMF),
        Some(Some(&InfoValue::Float(1.0)))
    );
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

#[test]
fn encodes_reference_calls_without_an_empty_ac_field() {
    let likelihood_schema =
        LikelihoodVcfSchema::new([(b"chr1".as_slice(), 5)], ["sample"]).unwrap();
    let likelihood = LikelihoodSite::new(
        0,
        0,
        Allele::new(&b"A"[..]).unwrap(),
        [Allele::new(&b"<*>"[..]).unwrap()],
        [1.0, 0.0],
        [SampleLikelihood::observed(
            Ploidy::new(2).unwrap(),
            [0, 3, 40],
            SampleEvidence::new(1, [1, 0], [40, 0]).unwrap(),
        )
        .unwrap()],
    )
    .unwrap();
    let called = crate::MultiallelicCaller::default()
        .call(&likelihood)
        .unwrap();
    let schema = CalledVcfSchema::from_likelihood(&likelihood_schema);
    let record = schema.encode_call(&called).unwrap();

    assert!(record.alternate_bases().as_ref().is_empty());
    assert!(record.info().get(AC).is_none());
    assert_eq!(record.info().get(AN), Some(Some(&InfoValue::Integer(2))));
    let mut data = Vec::new();
    let mut writer = vcf::io::Writer::new(&mut data);
    writer.write_header(schema.header()).unwrap();
    writer
        .write_variant_record(schema.header(), &record)
        .unwrap();
    let fields = std::str::from_utf8(&data)
        .unwrap()
        .lines()
        .last()
        .unwrap()
        .split('\t')
        .collect::<Vec<_>>();
    assert_eq!(fields[4], ".");
    assert_eq!(fields[7], "AN=2");
    assert_eq!(fields[9], "0/0:.:1:1:40:.:0");
}
