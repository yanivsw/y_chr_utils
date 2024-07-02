range_lower <- function(range) {
    vapply(range, function(x) {
        if (length(x) == 0)
            return(NA)
        as.numeric(x[1])
    }, numeric(1))
}

range_upper <- function(range) {
    vapply(range, function(x) {
        if (length(x) == 0)
            return(NA)
        if (length(x) == 1 && is.na(x))
            return(NA)
        as.numeric(x[2])
    }, numeric(1))
}

tree_metadata = data.frame(label = c(# Neandertals
                                     'Chagyrskaya2', 'Mezmaiskaya2, ~43 ka',
                                     # Africans
                                     'A00 (Shum Laka)', 'A00', 'A1b1b2b', 'A1b1a1a1', 'A1a', 'A0a1', 'B-M181',
                                     'B2b2', 'B2b1a1a', 'B2b1a1c', 'B2a', 'E1b1a1a1c2c3b', 'E2b1a1', 'E1a2a1b1',
                                     # non-Africans
                                     'C2b1c', 'D1a1a', 'G2b1', 'D1b2a2', 'L1a2', 'M1', 'N1', 'J2a1', 'J1a2b',
                                     'I2a2a1a2a2', 'O1b1a1a1a', 'O2a2b1a1a', 'Q1a2a1a1', 'R2a', 'R1a1a1b2',
                                     'K (Ust\'Ishim, ~44 ka)', 'P1 (Yana1, ~32 ka)', 'I2 (Loschbour, ~8 ka)'),
                           group = c(rep("Neandertal", 2), rep("African", 14), rep("non-African", 18)))

fix_tree_labels <- function(tr) {
    tr$tip.label[tr$tip.label=='HG02982'] = 'A0a1'
    tr$tip.label[tr$tip.label=='HG02666'] = 'A1a'
    tr$tip.label[tr$tip.label=='HGDP01406'] = 'A1b1b2b'
    tr$tip.label[tr$tip.label=='HGDP01029'] = 'A1b1a1a1'
    tr$tip.label[tr$tip.label=='HGDP00931'] = 'B-M181'
    tr$tip.label[tr$tip.label=='HGDP00478'] = 'B2a'
    tr$tip.label[tr$tip.label=='HGDP00453'] = 'B2b2'
    tr$tip.label[tr$tip.label=='HGDP00984'] = 'B2b1a1a'
    tr$tip.label[tr$tip.label=='HGDP00475'] = 'B2b1a1c'
    tr$tip.label[tr$tip.label=='HGDP01200'] = 'E1a2a1b1'
    tr$tip.label[tr$tip.label=='HGDP00908'] = 'E1b1a1a1c2c3b'
    tr$tip.label[tr$tip.label=='HGDP01031'] = 'E2b1a1'

    tr$tip.label[tr$tip.label=='HGDP00103'] = 'C2b1c'
    tr$tip.label[tr$tip.label=='HGDP01214'] = 'D1a1a'
    tr$tip.label[tr$tip.label=='HGDP00757'] = 'D1b2a2'
    tr$tip.label[tr$tip.label=='HGDP00213'] = 'G2b1'
    tr$tip.label[tr$tip.label=='HGDP00056'] = 'J2a1'
    tr$tip.label[tr$tip.label=='HGDP00057'] = 'J1a2b'
    tr$tip.label[tr$tip.label=='HGDP00127'] = 'I2a2a1a2a2'
    tr$tip.label[tr$tip.label=='HGDP01190'] = 'O1b1a1a1a'
    tr$tip.label[tr$tip.label=='HGDP01225'] = 'O2a2b1a1a'
    tr$tip.label[tr$tip.label=='HGDP01009'] = 'Q1a2a1a1'
    tr$tip.label[tr$tip.label=='HGDP00218'] = 'R2a'
    tr$tip.label[tr$tip.label=='HGDP00136'] = 'R1a1a1b2'
    tr$tip.label[tr$tip.label=='HGDP01298'] = 'L1a2'
    tr$tip.label[tr$tip.label=='HGDP01192'] = 'N1'
    tr$tip.label[tr$tip.label=='HGDP00549'] = 'M1'

    tr$tip.label[tr$tip.label=='Ust_Ishim']    = 'K (Ust\'Ishim, ~44 ka)'
    tr$tip.label[tr$tip.label=='Loschbour']    = 'I2 (Loschbour, ~8 ka)'
    tr$tip.label[tr$tip.label=='Yana1']        = 'P1 (Yana1, ~32 ka)'
    # tr$tip.label[tr$tip.label=='A00_I10871']   = 'A00 (Shum Laka, ~8 ka)'
    tr$tip.label[tr$tip.label=='A00_I10871']   = 'A00 (Shum Laka)'
    tr$tip.label[tr$tip.label=="Mezmaiskaya2"] = 'Mezmaiskaya2, ~43 ka'

    return(tr)
}