/*
 * This repository / codebase is Open Source and free for use and rewrite.
 */
package io.onclave.nsga.ii.datastructure;

import io.onclave.nsga.ii.api.Service;
import java.util.ArrayList;
import java.util.List;

/*这是一个生物染色体的模拟，它包含一个遗传密码，一个适应度，一个统治等级，一个唯一的ID，和一个被控制的染色体的列表。*/

/**
 * this is a simulation of a biological chromosome that contains a genetic code, a fitness value,
 * a domination rank, an unique ID, and a list of dominated chromosomes.
 * 
 * @author  Debabrata Acharya <debabrata.acharya@icloud.com>
 * @version 1.1
 * @since   0.1
 */
public class Chromosome {
    
    public Chromosome() {
        this(-Double.MIN_VALUE);
    }
    
    public Chromosome(final double fitness) {
        this(null, fitness, "");
    }
    
    public Chromosome(final Allele[] geneticCode, final double fitness, final String extraInfo) {
        this(geneticCode, fitness, extraInfo, 0);
        
    }
    
    public Chromosome(final Allele[] geneticCode, final double fitness, final String extraInfo, final int rank) {
        
        this.geneticCode = geneticCode;
        this.fitness = fitness;
        this.extraInfo = extraInfo;
        this.dominationRank = rank;
        this.uniqueID = Long.toString(System.currentTimeMillis()) + "-" + Integer.toString(Service.generateRandomInt());
        this.dominatedChromosomes = new ArrayList<>();
    }
    
    private Allele[] geneticCode;
    private double fitness;
    private String extraInfo;
    private int dominationRank = 0;
    private String uniqueID;
    private List<Chromosome> dominatedChromosomes;//被支配的染色体的列表

    public List<Chromosome> getDominatedChromosomes() {
        return dominatedChromosomes;
    }

    public void setDominatedChromosomes(List<Chromosome> dominatedChromosomes) {
        this.dominatedChromosomes = dominatedChromosomes;
    }

    public String getUniqueID() {//->41行
        return uniqueID;
    }

    public void setUniqueID(String uniqueID) {
        this.uniqueID = uniqueID;
    }

    public int getDominationRank() {
        return dominationRank;
    }

    public void setDominationRank(int dominationRank) {
        this.dominationRank = dominationRank;
    }

    public String getExtraInfo() {
        return extraInfo;
    }

    public void setExtraInfo(String extraInfo) {
        this.extraInfo = extraInfo;
    }

    public Allele[] getGeneticCode() {
        return geneticCode;
    }

    /**
     * the new fitness value is set as soon as a new genetic code is set for a chromosome.
     * @param   geneticCode     the genetic code that the chromosome carries.
     */
    public void setGeneticCode(Allele[] geneticCode) {
        this.geneticCode = geneticCode;
        this.setFitness(Service.calculateFitness(geneticCode));//！！！->Service552
    }

    public double getFitness() {
        return fitness;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }
    
    @Override
    public String toString() {
    	String geneticString = "";
    	for(Allele allele : geneticCode) geneticString += allele.getGene() ? "1" : "0";
    	return geneticString;
    }
}
