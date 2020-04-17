/*
 * This repository / codebase is Open Source and free for use and rewrite.
 */
package io.onclave.nsga.ii.api;
//合成随机的染色体群,Selection,Crossover,Mutation
import io.onclave.nsga.ii.configuration.Configuration;
import io.onclave.nsga.ii.datastructure.Allele;
import io.onclave.nsga.ii.datastructure.Chromosome;
import io.onclave.nsga.ii.datastructure.Population;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

//这是一个合成类，它执行许多底层工作（生物模拟），这些工作是由业务/控制器层的其他类抽象或封装的。

/**
 * This is the synthesis class that does many of the under-the-hood work (biological simulation) that is abstracted/encapsulated
 * by other classes at the business/controller layer.
 * 
 * @author  Debabrata Acharya <debabrata.acharya@icloud.com>
 * @version 1.1
 * @since   0.2
 */
public class Synthesis {
    
    private static final Random LOCAL_RANDOM = new Random();
    //这种方法合成了随机的染色体群，并为每条染色体提供了伪随机生成的遗传密码。
    /**
     * depending on the settings available in the Configuration.java file, this method synthesizes a
     * random population of chromosomes with pseudo-randomly generated genetic code for each chromosome.
     * 
     * @return  a randomly generated population
     */
    public static Population syntesizePopulation() {
        
        List<Chromosome> populace = new ArrayList<>();
        
        /**
         * the number of chromosomes in the population is received from the Configuration.java file
         */
        for(int i = 0; i < Configuration.getPOPULATION_SIZE(); i++) {
        	//底层synthesizeGeneticCode(final int length)已修改，产生的染色体已满足要求
        	Chromosome chromosome = new Chromosome();
        	do {
            chromosome.setGeneticCode(synthesizeGeneticCode(Configuration.getCHROMOSOME_LENGTH()));//set fitness
        	}while(populace.contains(chromosome));//保证染色体无重复
        	populace.add(chromosome);//GeneticCode_sizw=30bit
        }
        
        return new Population(populace);
    }
    
    /**
     * a child population of the same size as the parent is synthesized from the parent population
     * 
     * @param   parent  the parent population object
     * @return          a child population synthesized from the parent population
     */
    public static Population synthesizeChild(Population parent) {
    	//crossover,mutation已修改，产生满足要求的染色体
        Population child = new Population();
        List<Chromosome> populace = new ArrayList<>();
        
        /**
         * child chromosomes undergo crossover and mutation.
         * the child chromosomes are selected using binary tournament selection.
         * crossover returns an array of exactly two child chromosomes synthesized from two parent
         * chromosomes.
         */

        while(populace.size() < Configuration.getPOPULATION_SIZE()) {
        	Chromosome[] crossChromosome=new Chromosome[2];
        	do {//保证染色体不重
        		crossChromosome=crossover(binaryTournamentSelection(parent), binaryTournamentSelection(parent));
        	}while(populace.contains(crossChromosome[0])||populace.contains(crossChromosome[1]));
        	for(Chromosome chromosome : crossChromosome)
                populace.add(mutation(chromosome));
        }      	
       
        child.setPopulace(populace);
        
        return child;
    }
    
    /**
     * this is an implementation of basic binary tournament selection.
     * for a tournament of size t, select t individuals (randomly) from population and determine winner of
     * tournament with the highest fitness value.
     * in case of binary tournament selection, t = 2.
     * 
     * refer [https://stackoverflow.com/questions/36989783/binary-tournament-selection] for more information.
     * 
     * @param   population  the population from which a child chromosome is to be selected
     * @return              the selected child chromosome
     */
    private static Chromosome binaryTournamentSelection(Population population) {//锦标赛选择
        
              
    	Chromosome individual1 = population.getPopulace().get(LOCAL_RANDOM.nextInt(population.getPopulace().size()));                                      
        Chromosome individual2 = population.getPopulace().get(LOCAL_RANDOM.nextInt(population.getPopulace().size()));      
        
        if(individual1.getFitness() > individual2.getFitness()) return individual1; else return individual2;
    }
    
    /**
     * this is a basic implementation of uniform crossover where the crossover/break point is the middle
     * of the chromosomes. The genetic code of both the parent chromosomes are broken from the middle
     * and crossover is done to create two child chromosomes.
     * crossover probability is considered.
     * 
     * @param   chromosome1     the first parent chromosome taking part in crossover
     * @param   chromosome2     the second parent chromosome taking part in crossover
     * @return                  an array of exactly two child chromosomes synthesized from two parent chromosomes.
     */
    public static Chromosome[] crossover(Chromosome chromosome1, Chromosome chromosome2) {
        
        Allele[] geneticCode1 = new Allele[Configuration.getCHROMOSOME_LENGTH()];
        Allele[] geneticCode2 = new Allele[Configuration.getCHROMOSOME_LENGTH()];
        Allele[] chromosome1geneCode = chromosome1.getGeneticCode();
        Allele[] chromosome2geneCode = chromosome2.getGeneticCode();
        Chromosome[] childChromosomes = {new Chromosome(), new Chromosome()};
        //int breakPoint = Configuration.getCHROMOSOME_LENGTH() / 2;//交叉位置，单点交叉，染色体前一半不变，后一半交换
        
        Random random = new Random();	//均匀交叉
        int[] Point = new int[Configuration.getCHROMOSOME_LENGTH()];
       
        /**
         * generating a new random float value and if this value is less than equal to the
         * crossover probability mentioned in the Configuration file, then crossover occurs,
         * otherwise the parents themselves are copied as child chromosomes.
         */
        if(LOCAL_RANDOM.nextFloat() <= Configuration.getCROSSOVER_PROBABILITY()) {//达到概率交叉变异
        	
        		for(int i = 0;i < Configuration.getCHROMOSOME_LENGTH();i++)
            	Point[i] = random.nextInt(2);//产生[0,2)之间的随机整数
        	
        		for(int i = 0; i < Configuration.getCHROMOSOME_LENGTH(); i++) {
            
        			if(Point[i] == 0) {
        				geneticCode1[i] = chromosome1geneCode[i];
        				geneticCode2[i] = chromosome2geneCode[i];
        			} else {
        				geneticCode1[i] = chromosome2geneCode[i];
        				geneticCode2[i] = chromosome1geneCode[i];
        			}
        		}
    
            childChromosomes[0].setGeneticCode(geneticCode1);
            childChromosomes[1].setGeneticCode(geneticCode2);
           
        } else {//没有达到概率就不交叉
            childChromosomes[0] = chromosome1;
            childChromosomes[1] = chromosome2;
        }
        
        return childChromosomes;
    }
    
    /**
     * in this mutation operation implementation, a random bit-flip takes place.
     * a random float value is generated and if this value is less than equal to the mutation
     * probability defined in Configuration, then mutation takes place, otherwise the original
     * chromosome is returned.
     * 
     * @param   chromosome  the chromosome over which the mutation takes place
     * @return              the mutated chromosome
     */
    private static Chromosome mutation(Chromosome chromosome) {
        
        if(LOCAL_RANDOM.nextFloat() <= Configuration.getMUTATION_PROBABILITY()) {
            
            Allele[] geneticCode = chromosome.getGeneticCode();
            //do{
            geneticCode[LOCAL_RANDOM.nextInt(geneticCode.length)].bitFlip();//->Alleteȡ取反操作
            //}while(!valid(weights(geneticCode)));//变异附加条件
            
            chromosome.setGeneticCode(geneticCode);
        }
        
        return chromosome;
    }
    
    /**
     * a genetic code as an array of Alleles is synthesized.
     * refer Allele.java for more information.
     * 
     * @param   length  the required length of the genetic code.
     * @return          the synthesized genetic code.
     */
    public static Allele[] synthesizeGeneticCode(final int length) {//length=CHROMOSOME_LENGTH
        
    	Allele[] geneticCode = new Allele[length];
        
       //do{
        for(int i = 0; i < length; i++) geneticCode[i] = synthesizeAllele();
        //}while(!valid(weights(geneticCode)));//已修改，代码值必合格
        
        return geneticCode;
    }   
    
     //判断产生染色体是否满足要求
    /*public static boolean valid(double[] w) {//geneticCodeNumbers=10bit    final Allele[] geneticCode
      
    	double value1 = w[0],value2 = w[1];
    	//if(value1 + value2 <= 1.0) return true;
    	if(value1 + value2 == 1.0) return true;
        else return false;
    }*/
    
  //对染色体进行分解取值
    public static double[] weights(final Allele[] geneticCode) {
    	double[] weights = new double[2];
    	double value1 = 0,value2 = 0/*,value3 = 0*/;
        String binaryString = "";
        String binaryString1 ="" ;
        String binaryString2 = "";
        //String binaryString3 = "";
        int length = geneticCode.length/2;
        for(Allele bit : geneticCode) binaryString += bit.getGene() ? "1" : "0";
         
        
        binaryString1 = binaryString.substring(0, length);
        binaryString2 = binaryString.substring(length/*,2*length*/);
        //binaryString3 = binaryString.substring(2*length);
        
        for(int i = 0; i < binaryString1.length(); i++) if(binaryString1.charAt(i) == '1') value1 += Math.pow(2, binaryString1.length() - 1 - i);
        value1=value1/(Math.pow(2, binaryString1.length())-1);
        value1=(double)Math.round(value1*100)/100;
        
        for(int i = 0; i < binaryString2.length(); i++) if(binaryString2.charAt(i) == '1') value2 += Math.pow(2, binaryString2.length() - 1 - i);
        value2=value2/(Math.pow(2, binaryString2.length())-1);
        value2=(double)Math.round(value2*100)/100;
        
        
        /*for(int i = 0; i < binaryString3.length(); i++) if(binaryString3.charAt(i) == '1') value3 += Math.pow(2, binaryString3.length() - 1 - i);
        //二进制基因转换成是十进制
        value3=value3/(Math.pow(2, binaryString3.length()) - 1);//产生value=[0,1]->weightArray
        value3=(double)Math.round(value3*100)/100;//double后值保留两位->weightArray*/
        
        /*value3 = 1 - value1 - value2;
        value3=(double)Math.round( value3 * 100 ) / 100;*/
    	weights[0] = (double)Math.round((value1/(value1+value2))*100)/100;weights[1]=1-weights[0];
    	//weights[1] = value2/(value1+value2);//weights[2] = value3;
    	
    	return weights;
    }
    
    
    /**
     * an allele object with a randomly selected boolean gene value is synthesized.
     * 
     * @return  a randomly generated Allele object
     */
    public static Allele synthesizeAllele() {
        return new Allele(LOCAL_RANDOM.nextBoolean());
    }
    
    public static int[] Topsis(double[][] X) {
    	double[][] r = new double[Configuration.getPOPULATION_SIZE()][Configuration.getObjectives().size()];
    	double[][] v = new double[Configuration.getPOPULATION_SIZE()][Configuration.getObjectives().size()];;
    	double[] sum = {0,0},VL = {0,0},VF = {0,0};   	
    	double[] DL = new double[Configuration.getPOPULATION_SIZE()];
    	double[] DF = new double[Configuration.getPOPULATION_SIZE()];
    	
    	//��ƽ����
    	for(int j = 0;j <  Configuration.getObjectives().size();j++)
    		for(int i = 0;i < Configuration.getPOPULATION_SIZE();i++) {
    		sum[j] += Math.pow(X[i][j], 2);
    	}
    	
    	//求平方和
    	for(int i = 0;i < Configuration.getPOPULATION_SIZE();i++)
    		for(int j = 0;j < Configuration.getObjectives().size();j++) {
    			r[i][j] = X[i][j] / Math.sqrt(sum[j]);
    		}
    	
    	//标准化决策矩阵(ij)
    	for(int i = 0;i < Configuration.getPOPULATION_SIZE();i++)
    		for(int j = 0;j < Configuration.getObjectives().size();j++) {
    			v[i][j] = r[i][j]*X[i][j];
    		}
    	
    	//确定理想解和负理想解VL/VF
    	for(int j = 0;j <  Configuration.getObjectives().size();j++)
    		for(int i = 0;i < Configuration.getPOPULATION_SIZE();i++) {
    		if(VL[j] <= v[i][j]) VL[j] = v[i][j];
    		if(VF[j] >= v[i][j]) VF[j] = v[i][j];
    	}
    	
    	//计算每个备选方案到理想解与负理想解的距离DL(i),DF(i)
    	for(int i = 0;i < Configuration.getPOPULATION_SIZE();i++) 
    		for(int j = 0;j <  Configuration.getObjectives().size();j++) {
    			double sL = 0,sF = 0;
    			sL+= Math.pow((v[i][j] - VL[j]), 2); 
        		DL[i] = Math.sqrt(sL);
        		sF+= Math.pow((v[i][j] - VF[j]), 2); 
        		DF[i] = Math.sqrt(sF);
    		}
    	
    	//计算最接近理想解的方案C(i)
    	double[] C = new double[Configuration.getPOPULATION_SIZE()];
    	int[] Ci = new int[Configuration.getPOPULATION_SIZE()];
    	for(int i = 0;i < Configuration.getPOPULATION_SIZE();i++) {
    		C[i] = DF[i] / (DL[i] + DF[i]);
    		Ci[i] = i;
    	}
    	//对Ci排序
    	for(int i = 0;i < Configuration.getPOPULATION_SIZE()-1;i++)
    		for(int j = 1;j < Configuration.getPOPULATION_SIZE()-i-1;j++)
    			if(C[i] < C[i+j]) {
    				double temp;temp = C[i];C[i] = C[i+j];C[i+j] = temp;
    				int tempi;tempi = Ci[i];Ci[i] = Ci[i+j];Ci[i+j] = tempi;
    			}
    	
    	/*for(int i = 0;i < 10;i++) {
    		System.out.println(C[i] + " " + Ci[i]);
    		
    	}*/
    	return Ci;
    }
    
    public static void main(String[] args) {
    	
    	/*for(int j=0;j<100;j++) {
    	double w[]=weights(synthesizeGeneticCode(12));
    	System.out.println(w[0]+" "+w[1]+"\t=\t"+(w[0]+w[1])+"\n");}*/
    	
    	Population testpopulace = Service.nonDominatedPopulationSort(Synthesis.syntesizePopulation());
    	
    	for(int j=0;j<100;j++) {
    		double w[]=weights(binaryTournamentSelection(testpopulace).getGeneticCode());
    		System.out.println( w[0] + "\t" +w[1] +"\t=\t"+(w[0]+w[1])+"\n");
    		}
    	
    	
    	/*for(int j=0;j<50;j++) {
    		for(Chromosome chromosome : crossover(binaryTournamentSelection(testpopulace), binaryTournamentSelection(testpopulace))) {
    			double w[]=weights(chromosome.getGeneticCode());
    			System.out.println( w[0] + "\t" +w[1] +"\t=\t"+(w[0]+w[1])+"\n");
    		}
    	}*/
    	
    	/*for(int j=0;j<50;j++) {
    		double w[]=weights(mutation(testpopulace.getPopulace().get(j)).getGeneticCode());
    		System.out.println( w[0] + "\t" +w[1] +"\t=\t"+(w[0]+w[1])+"\n");
    	}*/

    } 
    
}
