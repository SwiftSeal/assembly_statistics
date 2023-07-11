mod utils;

use bio::io::fasta;
use std::io;
use crate::utils::*;

fn main() {
    let reader: fasta::Reader<io::BufReader<io::Stdin>> = fasta::Reader::new(io::stdin());

    let mut scaffolds: Vec<Scaffold> = Vec::new();

    for result in reader.records() {
        let record: fasta::Record = result.expect("Error during fasta record parsing");

        let scaffold: Scaffold = parse_scaffold(record.seq());

        scaffolds.push(scaffold);
    }

    let summary: Summary = summary_statistics(&scaffolds);

    print_summary(&summary);

}
