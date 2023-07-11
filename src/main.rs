struct Scaffold {
    a: u64,
    t: u64,
    g: u64,
    c: u64,
    length: u64,
    contig_lengths: Vec<u64>,
    contig_gcs: Vec<f64>,
}

struct Summary {
    n50: u64,
    l50: u64,
    n90: u64,
    l90: u64,
    scaffold_length: u64,
    scaffold_count: u64,
    contig_length: u64,
    contig_count: u64,
    largest_scaffold: u64,
    largest_contig: u64,
    sorted_contig_lengths: Vec<u64>,
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

impl std::fmt::Display for Summary {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "N50: {}\nL50: {}\nN90: {}\nL90: {}\nScaffold length: {}\nScaffold count: {}\nContig length: {}\nContig count: {}\nLargest scaffold: {}\nLargest contig: {}",
            self.n50,
            self.l50,
            self.n90,
            self.l90,
            self.scaffold_length,
            self.scaffold_count,
            self.contig_length,
            self.contig_count,
            self.largest_scaffold,
            self.largest_contig
        )
    }
}

fn unit_scaling(value: u64, decimal: usize) -> String {
    // This function takes the number of basepairs and returns a string with an appropriate unit. Units should be < 1000 and rounded to 1 decimal place
    if value < 1000 {
        format!("{} bp", value)
    } else if value < 1000000 {
        format!("{:.1$} kbp", value as f64 / 1000.0, decimal)
    } else if value < 1000000000 {
        format!("{:.1$} Mbp", value as f64 / 1000000.0, decimal)
    } else {
        format!("{:.1$} Gbp", value as f64 / 1000000000.0, decimal)
    }
}

fn parse_scaffold(sequence: &[u8]) -> Scaffold {
    let mut a: u64 = 0;
    let mut t: u64 = 0;
    let mut g: u64 = 0;
    let mut c: u64 = 0;
    let mut length: u64 = 0;
    let mut contig_lengths: Vec<u64> = Vec::new();
    let mut contig_gcs: Vec<f64> = Vec::new();

    let mut contig_length: u64 = 0;
    let mut contig_gc: i32 = 0;

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

fn summary_statistics(scaffolds: &[Scaffold]) -> Summary {
    // Scaffold statistics
    let scaffold_length: u64 = scaffolds.iter().map(|scaffold| scaffold.length).sum();
    let scaffold_count: u64 = scaffolds.len() as u64;
    let largest_scaffold: u64 = scaffolds
        .iter()
        .map(|scaffold: &Scaffold| scaffold.length)
        .max()
        .unwrap_or(0);

    // get sorted list of all contig lengths
    let mut contig_lengths: Vec<u64> = scaffolds
        .iter()
        .flat_map(|scaffold: &Scaffold| scaffold.contig_lengths.clone())
        .collect();

    contig_lengths.sort();

    let largest_contig: u64 = contig_lengths.last().copied().unwrap_or(0);

    let contig_length: u64 = contig_lengths.iter().sum();
    let contig_count: u64 = contig_lengths.len() as u64;
    let n50_threshold: u64 = (contig_length as f64 * 0.5) as u64;
    let n90_threshold: u64 = (contig_length as f64 * 0.9) as u64;

    // iterate over contig lengths in descending order and calculate N/L50 and N/L90. Once N50 has been calculated it will not be recalculated
    let mut accumulated_length: u64 = 0;
    let mut number_of_contigs: u64 = 0;

    let mut n50: u64 = 0;
    let mut l50: u64 = 0;
    let mut n90: u64 = 0;
    let mut l90: u64 = 0;
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

    Summary {
        n50,
        l50,
        n90,
        l90,
        scaffold_length,
        scaffold_count,
        contig_length,
        contig_count,
        largest_scaffold,
        largest_contig,
        sorted_contig_lengths: contig_lengths,
    }
}

fn print_summary(summary: &Summary) {
    println!("Scaffold count: {}", summary.scaffold_count);
    println!("Scaffold length: {}", unit_scaling(summary.scaffold_length, 1));
    println!("Largest scaffold: {}", unit_scaling(summary.largest_scaffold, 1));
    println!("Contig count: {}", summary.contig_count);
    println!("Contig length: {}", unit_scaling(summary.contig_length, 1));
    println!("Largest contig: {}", unit_scaling(summary.largest_contig, 1));
    println!("L/N50: {}/{}", summary.l50, unit_scaling(summary.n50, 1));
    println!("L/N90: {}/{}", summary.l90, unit_scaling(summary.n90, 1)); 

    let mut minimum_length: u64 = 1000;
    loop {
        let mut n_contigs: u64 = 0;
        let mut accumulated_length: u64 = 0;

        for length in summary.sorted_contig_lengths.iter() {
            if *length >= minimum_length {
                n_contigs += 1;
                accumulated_length += length;
            }
        }

        if n_contigs == 0 {
            break;
        }

        println!("{}\t{}\t{}", unit_scaling(minimum_length, 0), n_contigs, unit_scaling(accumulated_length, 1));

        minimum_length *= 10;
    }
}

fn main() {
    use bio::io::fasta;
    use std::io;
    use std::fs::File;
    use std::io::Write;

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
