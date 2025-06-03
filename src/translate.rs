/// Attempt to translate a DNA sequence to an amino acid sequence. If the sequence
/// length is not divisible by 3, the None variant is returned.
pub fn translate(seq: &[u8]) -> Option<String> {
    if seq.len() % 3 != 0 {
        return None;
    }

    let mut peptide = String::with_capacity(seq.len() / 3);

    'outer: for triplet in seq.chunks_exact(3) {
        for c in triplet {
            if !c.is_ascii() {
                peptide.push('X');
                continue 'outer;
            }
        }

        let c1 = ASCII_TO_INDEX[triplet[0] as usize];
        let c2 = ASCII_TO_INDEX[triplet[1] as usize];
        let c3 = ASCII_TO_INDEX[triplet[2] as usize];

        let amino_acid = if c1 == 4 || c2 == 4 || c3 == 4 {
            'X'
        } else {
            AA_TABLE_CANONICAL[c1][c2][c3]
        };

        peptide.push(amino_acid);
    }
    Some(peptide)
}

/// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
/// U is equivalent to T here
///
/// The 1st index picks the 4x4 block
/// The 2nd index picks the row
/// the 3rd index picks the column
static AA_TABLE_CANONICAL: [[[char; 4]; 4]; 4] = [
    [
        ['K', 'N', 'K', 'N'], // AAA, AAC, AAG, AAU/AAT
        ['T', 'T', 'T', 'T'], // ACA, ACC, ACG, ACU/ACT
        ['R', 'S', 'R', 'S'], // AGA, AGC, AGG, AGU/AGT
        ['I', 'I', 'M', 'I'], // AUA/ATA, AUC/ATC, AUG/ATG, AUU/ATT
    ],
    [
        ['Q', 'H', 'Q', 'H'], // CAA, CAC, CAG, CAU/CAT
        ['P', 'P', 'P', 'P'], // CCA, CCC, CCG, CCU/CCT
        ['R', 'R', 'R', 'R'], // CGA, CGC, CGG, CGU/CGT
        ['L', 'L', 'L', 'L'], // CUA/CTA, CUC/CTC, CUG/CTG, CUU/CTT
    ],
    [
        ['E', 'D', 'E', 'D'], // GAA, GAC, GAG, GAU/GAT
        ['A', 'A', 'A', 'A'], // GCA, GCC, GCG, GCU/GCT
        ['G', 'G', 'G', 'G'], // GGA, GGC, GGG, GGU/GGT
        ['V', 'V', 'V', 'V'], // GUA/GTA, GUC/GTC, GUG/GTG, GUU/GTT
    ],
    [
        ['*', 'Y', '*', 'Y'], // UAA/TAA, UAC/TAC, UAG/TAG, UAU/TAT
        ['S', 'S', 'S', 'S'], // UCA/TCA, UCC/TCC, UCG/TCG, UCU/TCT
        ['*', 'C', 'W', 'C'], // UGA/TGA, UGC/TGC, UGG/TGG, UGU/TGT
        ['L', 'F', 'L', 'F'], // UUA/TTA, UUC/TTC, UUG/TTG, UUU/TTT
    ],
];

/// Maps an ASCII character to array index
///
/// A = 65, a = 97  => 0
/// C = 67, c = 99  => 1
/// G = 71, g = 103 => 2
/// T = 84, t = 116 => 3
/// U = 85, u = 117 => 3
static ASCII_TO_INDEX: [usize; 128] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 48-63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 64-79    (65 = A, 67 = C, 71 = G)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 80-95    (84 = T, 85 = U)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 96-111   (97 = a, 99 = c, 103 = g)
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 112-127  (116 = t, 117 = u)
];
