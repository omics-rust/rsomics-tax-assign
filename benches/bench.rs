use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_tax_assign(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-tax-assign");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let reads = manifest.join("tests/golden/reads.fq");
    let db = manifest.join("tests/golden/empty_db.tsv");
    c.bench_function("rsomics-tax-assign golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([reads.to_str().unwrap(), "-d", db.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_tax_assign);
criterion_main!(benches);
