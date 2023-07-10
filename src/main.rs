struct Scaffold {
    a: u64,
    t: u64,
    g: u64,
    c: u64,
    length: u64,
    contig_lengths: Vec<u64>,
    contig_gcs: Vec<f64>,
}

impl std::fmt::Display for Scaffold {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.a,
            self.t,
            self.g,
            self.c,
            self.length,
            self.contig_lengths
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(","),
            self.contig_gcs
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(",")
        )
    }
}

fn parse_scaffold(sequence: &[u8]) -> Scaffold {
    let mut a = 0;
    let mut t = 0;
    let mut g = 0;
    let mut c = 0;
    let mut length = 0;
    let mut contig_lengths = Vec::new();
    let mut contig_gcs = Vec::new();

    let mut contig_length = 0;
    let mut contig_gc = 0;

    for base in sequence {
        match base {
            b'A' => {
                a += 1;
                contig_length += 1;
            }
            b'T' => {
                t += 1;
                contig_length += 1;
            }
            b'G' => {
                g += 1;
                contig_length += 1;
                contig_gc += 1;
            }
            b'C' => {
                c += 1;
                contig_length += 1;
                contig_gc += 1;
            }
            b'N' => {
                if contig_length > 0 {
                    contig_lengths.push(contig_length);
                    contig_gcs.push(contig_gc as f64 / contig_length as f64);
                }
                contig_length = 0;
                contig_gc = 0;
            }
            _ => (),
        }
        length += 1;
    }

    if contig_length > 0 {
        contig_lengths.push(contig_length);
        contig_gcs.push(contig_gc as f64 / contig_length as f64);
    }

    Scaffold {
        a,
        t,
        g,
        c,
        length,
        contig_lengths,
        contig_gcs,
    }
}

fn calculate_ln50(scaffolds: &[Scaffold]) -> u64 {
    // Calculate the N50 and L50 of the scaffolds
    let mut lengths: Vec<u64> = scaffolds
        .iter()
        .map(|scaffold| scaffold.length)
        .collect();
    lengths.sort();

    let mut total_length = 0;
    for length in lengths.iter().rev() {
        total_length += length;
        if total_length >= scaffolds.iter().map(|scaffold| scaffold.length).sum::<u64>() / 2 {
            return *length;
        }
    }
    0
}

fn summary_statistics(scaffolds: &[Scaffold]) {
    let mut n50 = 0;
    let mut l50 = 0;
    let mut n90 = 0;
    let mut l90 = 0;
    let mut scaffold_length = 0;
    let mut scaffold_count = 0;
    let mut contig_length = 0;
    let mut contig_count = 0;
    let mut largest_scaffold = 0;
    let mut largest_contig = 0;

    // get sorted list of all contig lengths
    let mut contig_lengths: Vec<u64> = scaffolds
        .iter()
        .flat_map(|scaffold| scaffold.contig_lengths.clone())
        .collect();
    
    contig_lengths.sort();

    let mut accumulated_length = 0;
    let mut number_of_contigs = 0;

    let total_length: u64 = contig_lengths.iter().sum();
    let n50_threshold = (total_length as f64 * 0.5) as u64;
    let n90_threshold = (total_length as f64 * 0.9) as u64;

    // iterate over contig lengths in descending order and calculate N/L50 and N/L90. Once N50 has been calculated it will not be recalculated
    for length in contig_lengths.iter().rev() {
        accumulated_length += length;
        number_of_contigs += 1;

        if n50 == 0 && accumulated_length >= n50_threshold {
            n50 = *length;
            l50 = number_of_contigs;
        }

        if n90 == 0 && accumulated_length >= n90_threshold {
            n90 = *length;
            l90 = number_of_contigs;
            break;
        }
    }
        
}


fn main() {
    use bio::io::fasta;
    use std::io;

    let reader = fasta::Reader::new(io::stdin());

    let mut scaffolds: Vec<Scaffold> = Vec::new();

    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");

        let scaffold = parse_scaffold(record.seq());

        scaffolds.push(scaffold);
    }

    summary_statistics(&scaffolds);
}
