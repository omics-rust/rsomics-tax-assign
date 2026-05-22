use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-tax-assign"))
}
fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

/// Canonical k-mer hashes of `seq` via the same engine the tool uses.
fn kmers(seq: &[u8], k: usize) -> Vec<u64> {
    rsomics_kmer::KmerIter::new(seq, k, true)
        .unwrap()
        .flatten()
        .collect()
}

#[test]
fn output_has_correct_format() {
    let out = Command::new(ours())
        .arg(golden("reads.fq"))
        .args(["-d", &golden("empty_db.tsv"), "-k", "5"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    for line in s.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 4, "each line should have 4 fields");
        assert!(parts[0] == "C" || parts[0] == "U", "first field C or U");
    }
}

// Correctness: a read whose k-mers populate the DB under a taxon must be
// classified to that taxon; an unrelated read must be unclassified. The DB is
// built with the same k-mer engine the tool classifies with (round-trip).
#[test]
fn classifies_reads_to_source_taxon() {
    let k = 7;
    let refseq: &[u8] = b"ACGTACAGTCGATCGGATCCAGTACGTAGCATGCATCGATCGTTAGCC";
    let unrelated: &[u8] = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";

    let dir = std::env::temp_dir().join("rsomics-tax-assign-compat");
    let _ = std::fs::create_dir_all(&dir);
    let dbp = dir.join("db.tsv");
    let mut db = String::new();
    for km in kmers(refseq, k) {
        db.push_str(&format!("{km}\t42\n"));
    }
    std::fs::write(&dbp, db).unwrap();

    let rs = std::str::from_utf8(refseq).unwrap();
    let us = std::str::from_utf8(unrelated).unwrap();
    let reads = format!(
        "@r1\n{rs}\n+\n{}\n@r2\n{us}\n+\n{}\n",
        "I".repeat(rs.len()),
        "I".repeat(us.len())
    );
    let readsp = dir.join("reads.fq");
    std::fs::write(&readsp, reads).unwrap();

    let out = Command::new(ours())
        .arg(&readsp)
        .args(["-d"])
        .arg(&dbp)
        .args(["-k", "7"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = s.lines().collect();

    let r1 = lines.iter().find(|l| l.contains("r1")).unwrap();
    assert!(
        r1.starts_with("C\tr1\t42\t"),
        "read from the reference must classify to taxon 42: {r1}"
    );
    let r2 = lines.iter().find(|l| l.contains("r2")).unwrap();
    assert!(
        r2.starts_with("U\tr2\t"),
        "unrelated read must be unclassified: {r2}"
    );
}
