pub struct Scaffold {
    a: u64,
    t: u64,
    g: u64,
    c: u64,
    length: u64,
    contig_lengths: Vec<u64>,
    contig_gcs: Vec<f64>,
}

pub struct Summary {
    gc: f64,
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

pub fn unit_scaling(value: u64, decimal: usize) -> String {
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

pub fn parse_scaffold(sequence: &[u8]) -> Scaffold {
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

pub fn summary_statistics(scaffolds: &[Scaffold]) -> Summary {
    // Scaffold statistics
    let scaffold_length: u64 = scaffolds.iter().map(|scaffold| scaffold.length).sum();
    let scaffold_count: u64 = scaffolds.len() as u64;
    let largest_scaffold: u64 = scaffolds
        .iter()
        .map(|scaffold: &Scaffold| scaffold.length)
        .max()
        .unwrap_or(0);

    // calculate mean GC content
    let contig_gcs: Vec<f64> = scaffolds.iter().flat_map(|scaffold: &Scaffold| scaffold.contig_gcs.clone()).collect();
    let gc: f64 = contig_gcs.iter().sum::<f64>() / contig_gcs.len() as f64;

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
        gc,
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

pub fn print_summary(summary: &Summary) {
    use termion::style;
    println!("   {bold}Assembly size:{reset} {scaffold_length}", bold = style::Bold, reset = style::Reset, scaffold_length = unit_scaling(summary.scaffold_length, 1));
    println!(" {bold}Mean GC content:{reset} {gc:.1}%", bold = style::Bold, reset = style::Reset, gc = summary.gc * 100.0);
    println!("  {bold}Scaffold count:{reset} {scaffold_count}", bold = style::Bold, reset = style::Reset, scaffold_count = summary.scaffold_count);
    println!("  {bold}Scaffold L/N50:{reset} {l50}/{n50}", bold = style::Bold, reset = style::Reset, l50 = summary.l50, n50 = unit_scaling(summary.n50, 1));
    println!("  {bold}Scaffold L/N90:{reset} {l90}/{n90}", bold = style::Bold, reset = style::Reset, l90 = summary.l90, n90 = unit_scaling(summary.n90, 1));
    println!("{bold}Largest scaffold:{reset} {scaffold_length}", bold = style::Bold, reset = style::Reset, scaffold_length = unit_scaling(summary.largest_scaffold, 1));
    println!("    {bold}Contig count:{reset} {contig_count}", bold = style::Bold, reset = style::Reset, contig_count = summary.contig_count);
    println!("     {bold}Contig size:{reset} {contig_length}", bold = style::Bold, reset = style::Reset, contig_length = unit_scaling(summary.contig_length, 1));
    println!("  {bold}Largest contig:{reset} {largest_contig}", bold = style::Bold, reset = style::Reset, largest_contig = unit_scaling(summary.largest_contig, 1));

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

        println!(">{}\t{}\t{}", unit_scaling(minimum_length, 0), n_contigs, unit_scaling(accumulated_length, 1));

        minimum_length *= 10;
    }
}
