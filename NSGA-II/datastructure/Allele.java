/*
 * This repository / codebase is Open Source and free for use and rewrite.
 */
package io.onclave.nsga.ii.datastructure;

//这是对一个包含基因值的生物染色体中的等位基因的模拟。
//一组等位基为染色体创造了遗传密码。
/**
 * this is a simulation of an allele in a biological chromosome that contains a gene value.
 * an array of alleles create the genetic code for the chromosome.
 * 
 * @author  Debabrata Acharya <debabrata.acharya@icloud.com>
 * @version 1.0
 * @since   0.1
 */
public class Allele {
    
    public Allele() {
        this(false);
    }
    
    public Allele(final boolean gene) {
        this.gene = gene;
    }
    
    private boolean gene;

    public boolean getGene() {
        return gene;
    }

    public void setGene(boolean gene) {
        this.gene = gene;
    }
    
    public void bitFlip() {//比特反转
        this.gene = !this.gene;
    }
}
