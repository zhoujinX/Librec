/*
 * This repository / codebase is Open Source and free for use and rewrite.
 */
package io.onclave.nsga.ii.api;

import io.onclave.nsga.ii.Interface.IObjectiveFunction;
import io.onclave.nsga.ii.configuration.Configuration;
import io.onclave.nsga.ii.datastructure.Allele;
import io.onclave.nsga.ii.datastructure.Chromosome;
import io.onclave.nsga.ii.datastructure.ParetoObject;
import io.onclave.nsga.ii.datastructure.Population;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

/**
 * This is the service class that does most of the under-the-hood work that is abstracted/encapsulated
 * by other classes at the business/controller layer.
 * 
 * @author  Debabrata Acharya <debabrata.acharya@icloud.com>
 * @version 1.1
 * @since   0.1
 */
public class Service {
    
    /**
     * this is an implementation of the fast non-dominated sorting algorithm as defined in the
     * NSGA-II paper [DOI: 10.1109/4235.996017] Section III Part A.
     * 
     * @param   population  the population object that needs to undergo fast non-dominated sorting algorithm
     * @return  a HashMap with an integer key that labels the ranks and a list of chromosomes as values that clot chromosomes of same rank
     */
	static HashMap<Chromosome,Double> ObjFuction0 = new HashMap<>();//存储染色体对应的目标函数值
	static HashMap<Chromosome,Double> ObjFuction1 = new HashMap<>();
    public static HashMap<Integer, List<Chromosome>> fastNonDominatedSort(Population population) {
        //paretoFront作为最终返回，Integer为非支配等级，其对应为List<Chromosome>染色体集合
        HashMap<Integer, List<Chromosome>> paretoFront = new HashMap<>();//用来分类存储非支配集F1,F2,...Fn
        List<Chromosome> singularFront = new ArrayList<>();//临时存储个体chromosome.getDominationRank() == 0，即支配等级为0的个体
        List<Chromosome> populace = population.getPopulace();
//        System.out.println("计算目标函数");
        List<IObjectiveFunction> objectives = Configuration.getObjectives();
        for(Chromosome chromosome : populace) {//首先计算出染色体对应的目标函数值，保证每次迭代只需计算一次
        		ObjFuction0.put(chromosome, objectives.get(0).objectiveFunction(chromosome));
        		ObjFuction1.put(chromosome, objectives.get(1).objectiveFunction(chromosome));
        }
//        System.out.println("目标函数计算结束");
        /**
         * iterating over each chromosome of the population
         */
        for(Chromosome chromosome : populace) {
        	
            /**
             * an initial domination rank of 0 is set for each chromosome and a blank list is set for the number of
             * chromosomes that the present chromosome dominates.
             */
            chromosome.setDominationRank(0);//初始化操作，令所有染色体支配等级为0，并设置其支配染色体的集合
            chromosome.setDominatedChromosomes(new ArrayList<>());
            
            /**
             * for each chromosome, the program iterates over all the other remaining chromosomes to find which other
             * chromosomes are dominated by this chromosome and vice versa.
             */
            for (Chromosome competitor : populace) if(!competitor.equals(chromosome)) {
                
                /**
                 * if the present chromosome dominates the competitor, then:
                 *      i:   check if the competitor already exists in the list of dominated chromosomes of the present chromosome.
                 *     ii:   if the competitor does not exist within the list, then add it to the list of dominated chromosomes
                 *           of the present chromosome.
                 * else, if the competitor dominates the present chromosome, then increment the domination rank of the present
                 * chromosome by one.
                 */
                if(dominates(chromosome, competitor)) {//计算目标函数值
                    if(!chromosome.getDominatedChromosomes().contains(competitor)) chromosome.getDominatedChromosomes().add(competitor);
                } else if(dominates(competitor, chromosome)) chromosome.setDominationRank(chromosome.getDominationRank() + 1);
            }
            
            /**
             * if the domination rank of the present chromosome is 0, it means that this chromosome is a non-dominated chromosome
             * and hence it is added to the clot of chromosomes that are also non-dominated.
             */
            if(chromosome.getDominationRank() == 0&&!singularFront.contains(chromosome)) singularFront.add(chromosome);
        }
        
        /**
         * the first clot of non-dominated chromosomes is added to the HashMap with rank label 1.
         */
        paretoFront.put(1, singularFront);//将本次操作的支配等级为0的染色体集合（前沿面）加上非支配等级一起放到paretoFront中
        
        int i = 1;
        List<Chromosome> previousFront = paretoFront.get(i);//F1
        List<Chromosome> nextFront = new ArrayList<>();//F2,为previousFront中染色体所支配的染色体集合，作为下一个前沿面
        
        /**
         * the current/previous ranked clot of chromosomes with rank i is iterated over to find the next clot of chromosomes
         * with rank (i+1)
         */
        while(previousFront != null && !previousFront.isEmpty()) {
            
            /**
             * iterating over each chromosome from the previous clot of chromosomes ranked i.
             */
            for(Chromosome chromosome : previousFront) {
                
                /**
                 * iterating over each of the dominated chromosomes from the present chromosome of rank i.
                 */
                for(Chromosome recessive : chromosome.getDominatedChromosomes()) {
                    
                    /**
                     * if the domination rank of the current recessive chromosome in consideration is not 0, then
                     * decrement it's rank by 1.
                     * if the domination rank of the current recessive chromosome in consideration is 0, then add
                     * it to the next front [clot of chromosomes that belong to rank (i+1)].
                     */
                    if(recessive.getDominationRank() != 0) recessive.setDominationRank(recessive.getDominationRank() - 1);
                    if(recessive.getDominationRank() == 0) if(!nextFront.contains(recessive)) nextFront.add(recessive);
                }//F1所支配的染色体集合支配等级-1，由F1得到F2
            }
            
            /**
             * this code snippet ensures "rank jumps" to create all the possible rank lists from the parent
             * population.
             * new ranks are created only when there are recessive chromosomes with domination rank = 1 which are
             * decremented to domination rank 0 and then added to the next front.
             * but, due to the randomness of the algorithm, situation may occur such that even after decrementing all recessive
             * chromosome domination ranks by 1, none have domination rank 0 and hence the next front remains empty.
             * to ensure that all recessive chromosomes are added to some rank list, the program jumps domination ranks
             * of each recessive chromosome by decrementing domination rank by 1 until at least one of them reaches a
             * domination rank count of 0 and then that recessive chromosome is added to the next front.
             * if the next front is empty and the previous front has at least one dominated chromosome:
             *      i:  find the minimum rank among all the recessive chromosomes available:
             *              1:  iterate over all the chromosomes of the previous front
             *              2:  while the chromosomes have no dominated chromosomes with rank 0:
             *                      a:  iterate over all the recessive chromosomes of the current chromosome
             *                      b:  if the minimum rank is greater than the dominated rank of the present recessive,
             *                          mark this as the minimum rank recorded among all recessive chromosomes available.
             *              3:  end while
             *     ii:  iterate over all the chromosomes of the previous front
             *              1: while the chromosomes have no dominated chromosomes with rank 0:
             *                      a:  iterate over all the dominated chromosomes of the current chromosome
             *                      b:  if the domination rank of the recessive chromosome is not 0, then decrement the
             *                          domination count by value of minimum rank.
             *                      c:  if the domination rank is 0, then add it to the next front.
             *              2:  end while
             * 这里有一种情况，就是F1所支配的染色体集合中，每个染色体支配等级-1后不等于0，则无法加入F2中，此时F2就为空集
             */
            if(nextFront.isEmpty() && !isDominatedChromosomesEmpty(previousFront)) {//nextFront=F2
                
                int minimumRank = -1;
                
                for(Chromosome chromosome : previousFront)
                    while(hasRecessiveRankGreaterThanZero(chromosome))//检查输入染色体的任何一个支配染色体的支配秩是否为0，如果至少有一个支配染色体的支配秩为0，则返回true
                        for(Chromosome recessive : chromosome.getDominatedChromosomes())
                            if((minimumRank == -1) || minimumRank > recessive.getDominationRank()) minimumRank = recessive.getDominationRank();
                
                if(minimumRank != -1) for(Chromosome chromosome : previousFront)
                    while(hasRecessiveRankGreaterThanZero(chromosome)) for(Chromosome recessive : chromosome.getDominatedChromosomes()) {
                            if(recessive.getDominationRank() != 0) recessive.setDominationRank(recessive.getDominationRank() - minimumRank);
                            if(recessive.getDominationRank() == 0) if(!nextFront.contains(recessive)) nextFront.add(recessive);
                    }
            }
            
            /**
             * if the next front calculated is not empty, then it is added to the ranked HashMap data-structure
             * with the rank (i+1), else all chromosomes are sorted into some rank or the other and the program
             * breaks out of the loop.
             */
            if(!nextFront.isEmpty()) paretoFront.put(++i, nextFront); else break;
            
            /**
             * the next front (i) calculated is marked as the previous front for the next iteration (i+1) and
             * an empty next front is created.
             */
            previousFront = nextFront;
            nextFront = new ArrayList<>();
        }

        return paretoFront;//paretoFront(F1,F2,F2,...,Fn)
    }
    
    /**
     * this is the implementation of the crowding distance assignment algorithm as defined in the
     * NSGA-II paper [DOI: 10.1109/4235.996017] Section III Part B.
     * this ensures diversity preservation.
     * 
     * @param   singularFront   a list of chromosomes whose crowding distances are to be calculated
     * @return                  a list of ParetoObjects with assigned crowding distances. [Refer ParetoObject.java for more information]
     */
    public static List<ParetoObject> crowdingDistanceAssignment(List<Chromosome> singularFront) {
        //计算某一解集Fi的拥挤距离，选出前几个拥挤距离较大的解，保证多样性
        int i = 0;
        int end = singularFront.size() - 1;
        Double maxObjectiveValue;
        Double minObjectiveValue;
        List<IObjectiveFunction> objectives = Configuration.getObjectives();
        List<ParetoObject> singlePareto = new ArrayList<>();
        
        /**
         * for each chromosome in the input list, a new ParetoObject with an initial crowding distance of 0
         * is created and added to the list of ParetoObjects that are to be returned.
         */																//ParetoObject(chromosome, 0f)拥挤距离初始为float 0
        for(Chromosome chromosome : singularFront) singlePareto.add(i++, new ParetoObject(chromosome, 0f));
        
        /**
         * iterating over each of the objective functions set [refer Configuration.java for more information],
         * the ParetoObject list is sorted according to the objective functions and the first and last ParetoObjects
         * are set a crowding distance of infinity.
         */
        for(int j=0;j<2;j++) {
        	IObjectiveFunction objective=objectives.get(j);
            maxObjectiveValue = null;
            minObjectiveValue = null;
            singlePareto = sort(singlePareto, objective);
            
            //第一个染色体和最后一个染色体拥挤距离为无穷大
            singlePareto.get(0).setCrowdingDistance(Double.MAX_VALUE);
            singlePareto.get(end).setCrowdingDistance(Double.MAX_VALUE);
            
            /**
             * the max and min objective values are calculated according to the present objective function
             */
            for(ParetoObject paretoObject : singlePareto) {
                if(j==1) {
                	if((maxObjectiveValue == null) || (maxObjectiveValue < ObjFuction1.get(paretoObject.getChromosome()))) maxObjectiveValue = ObjFuction1.get(paretoObject.getChromosome());
                	if((minObjectiveValue == null) || (minObjectiveValue > ObjFuction1.get(paretoObject.getChromosome()))) minObjectiveValue = ObjFuction1.get(paretoObject.getChromosome());
                	}
                if(j==0) {
                    if((maxObjectiveValue == null) || (maxObjectiveValue < ObjFuction0.get(paretoObject.getChromosome()))) maxObjectiveValue = ObjFuction0.get(paretoObject.getChromosome());
                    if((minObjectiveValue == null) || (minObjectiveValue > ObjFuction0.get(paretoObject.getChromosome()))) minObjectiveValue = ObjFuction0.get(paretoObject.getChromosome());
                    }
                
            }
            
            /**
             * the crowding distance of all ParetoObjects are calculated and assigned except the first and last ParetoObjects
             * that have infinite crowding distance
             *///除去第一个和最后一个染色体的拥挤距离为无穷大，其他中间染色体计算拥挤距离
            for(i = 2; i < end; i++) singlePareto.get(i).setCrowdingDistance(calculateCrowdingDistance(singlePareto,
                                                                                                        i,
                                                                                                        j,
                                                                                                        objective,
                                                                                                        maxObjectiveValue,
                                                                                                        minObjectiveValue));
        }
        
        return singlePareto; 
   }
    
    /**
     * this method sorts a list of ParetoObjects based on the Crowd-Comparison Operator using the domination
     * rank and crowding distance as discussed in the NSGA-II paper [DOI: 10.1109/4235.996017] Section III Part B.
     * 
     * @param   singleFront     a list of ParetoObjects that are to be sorted according to their crowding distance
     * @return                  a list of sorted ParetoObjects
     */
    public static List<ParetoObject> crowdComparisonSort(List<ParetoObject> singleFront) {
        
        int index = -1;
        List<ParetoObject> sortedFront = new ArrayList<>();
        ParetoObject presentParetoObject;
        ParetoObject competitor;
        
        /**
         * all the ParetoObjects are, at first, marked as false for crowding distance sorted.
         */
        singleFront.stream().forEach((paretoObject) -> { paretoObject.setCrowdingDistanceSorted(false); });
        
        /**
         * iterating over each ParetoObject in the singular front input:
         *  i:  the i-th ParetoObject is marked as presentParetoObject
         * ii:  if the presentParetoObject is not already sorted by crowding distance:
         *          1:  iterate over the rest of the ParetoObjects in the input list as competitors that are
         *              not already sorted using crowding distance
         *          2:  compare the i-th and the j-th chromosome using the crowd comparison operator:
         *                  a: for different ranks, choose the one with the lower (better) rank.
         *                  b: for same rank, choose the one which has lower crowding distance.
         *          3:  if competitor dominates the i-th chromosome, then mark competitor as presentParetoObject
         *          4:  continue until i-th chromosome is compared to all competitors.
         *          5:  mark the presentParetoObject as already sorted by crowding distance
         *          6:  add presentParetoObject into list of sorted front with an incremented index
         */
        for(int i = 0; i < singleFront.size(); i++) {
            
            presentParetoObject = singleFront.get(i);
            
            if(!presentParetoObject.isCrowdingDistanceSorted()) {
                
                for(int j = 0; j < singleFront.size(); j++) {

                    competitor = singleFront.get(j);
                    
                    if(!competitor.isCrowdingDistanceSorted()) {
                        
                        double dominationRank = presentParetoObject.getChromosome().getDominationRank();
                        double competingDominationRank = competitor.getChromosome().getDominationRank();
                        double crowdingDistance = presentParetoObject.getCrowdingDistance();
                        double competingCrowdingDistance = competitor.getCrowdingDistance();

                        if(i != j) if((dominationRank > competingDominationRank) || ((dominationRank == competingDominationRank) && (crowdingDistance < competingCrowdingDistance))) presentParetoObject = competitor;
                    }
                }
                
                presentParetoObject.setCrowdingDistanceSorted(true);
                sortedFront.add(++index, presentParetoObject);
            }
        }
        
        return sortedFront;
    }
    
    /**
     * this method is not implemented, as it is not absolutely necessary for this algorithm to work.
     * is kept if implementation is needed in future.
     * returns the same unsorted parent population as of now.
     * 
     * @param   population  the population that is to be sorted
     * @return              a sorted population
     */
    public static Population nonDominatedPopulationSort(Population population) {
        
        //--TO-DO--
        
        return population;
    }
    
    /**
     * this method checks whether competitor1 dominates competitor2.
     * requires that none of the values of the objective functions using competitor1 is smaller
     * than the values of the objective functions using competitor2.
     * at least one of the values of the objective functions using competitor1 is greater than
     * the corresponding value of the objective functions using competitor2.
     * 
     * @param   competitor1     the chromosome that may dominate
     * @param   competitor2     the chromosome that may be dominated
     * @return                  boolean logic whether competitor1 dominates competitor2.
     */
     
    public static boolean dominates(final Chromosome competitor1, final Chromosome competitor2) {
    	
        /**
         * getting the list of configured objectives from Configuration.java
         */
    	if((ObjFuction0.get(competitor1)>ObjFuction0.get(competitor2))&&(ObjFuction1.get(competitor1)>ObjFuction1.get(competitor2)))
    		return true;
    	return false;
    }
    
    /**
     * the list is first converted to an array data-structure and then a randomized quick sort
     * algorithm is followed.
     * the resulting sorted array is again converted to a List data-structure before returning.
     * 
     * @param   singlePareto    the list of ParetoObjects that are to be sorted.
     * @param   objective       the objective function using which the ParetoObjects are sorted.
     * @return                  sorted list of ParetoObjects.
     */
    private static List<ParetoObject> sort(List<ParetoObject> singlePareto, IObjectiveFunction objective) {
        
        ParetoObject[] paretoArray = new ParetoObject[singlePareto.size()];
        singlePareto.toArray(paretoArray);
        
        randomizedQuickSort(paretoArray, 0, paretoArray.length - 1, objective);
        
        return (new ArrayList<>(Arrays.asList(paretoArray)));
    }
    
    /**
     * refer [https://jordanspencerwu.github.io/randomized-quick-sort/] for more details on randomized
     * quick sort algorithm.
     * 
     * @param   paretoArray     the array to be sorted
     * @param   head            the pointer/position of the head element
     * @param   tail            the pointer/position of the tail element
     * @param   objective       the objective function depending on which the sort is to take place
     * @return                  the pivot index.
     */
    private static int partition(ParetoObject[] paretoArray, int head, int tail, IObjectiveFunction objective) {
        
        ParetoObject pivot = paretoArray[tail];
        int i = head - 1;
        
        for(int j = head; j <= (tail - 1); j++) {
            
            if(objective.objectiveFunction(paretoArray[j]) <= objective.objectiveFunction(pivot)) {
                
                i++;
                ParetoObject temporary = paretoArray[i];
                paretoArray[i] = paretoArray[j];
                paretoArray[j] = temporary;
            }
        }
        
        ParetoObject temporary = paretoArray[i + 1];
        paretoArray[i + 1] = paretoArray[tail];
        paretoArray[tail] = temporary;
        
        return (i + 1);
    }
    
    /**
     * refer [https://jordanspencerwu.github.io/randomized-quick-sort/] for more details on randomized
     * quick sort algorithm.
     * 
     * @param   paretoArray     the array to be sorted
     * @param   head            the pointer/position of the head element
     * @param   tail            the pointer/position of the tail element
     * @param   objective       the objective function depending on which the sort is to take place
     * @return                  the random partition position index.
     */
    private static int randomizedPartition(ParetoObject[] paretoArray, int head, int tail, IObjectiveFunction objective) {
        
        int random = ThreadLocalRandom.current().nextInt(head, tail + 1);
        
        ParetoObject temporary = paretoArray[head];
        paretoArray[head] = paretoArray[random];
        paretoArray[random] = temporary;
        
        return partition(paretoArray, head, tail, objective);
    }
    
    /**
     * refer [https://jordanspencerwu.github.io/randomized-quick-sort/] for more details on randomized
     * quick sort algorithm.
     * 
     * @param   paretoArray     the array to be sorted
     * @param   head            the pointer/position of the head element
     * @param   tail            the pointer/position of the tail element
     * @param   objective       the objective function depending on which the sort is to take place
     */
    private static void randomizedQuickSort(ParetoObject[] paretoArray, int head, int tail, IObjectiveFunction objective) {
        
        if(tail < head) {
            
            int pivot = randomizedPartition(paretoArray, head, tail, objective);
            
            randomizedQuickSort(paretoArray, head, pivot - 1, objective);
            randomizedQuickSort(paretoArray, pivot + 1, tail, objective);
        }
    }
    
    /**
     * implementation of crowding distance calculation as defined in NSGA-II paper
     * [DOI: 10.1109/4235.996017] Section III Part B.
     * 
     * I[i]distance = I[i]distance + (I[i+1].m - I[i-1].m)/(f-max - f-min)
     * 
     * I[i]distance = crowding distance of the i-th individual
     * I[i+1].m = m-th objective function value of the (i+1)-th individual
     * I[i-1].m = m-th objective function value of the (i-1)-th individual
     * f-max, f-min = maximum and minimum values of the m-th objective function
     * 
     * @param   singlePareto            the list of ParetoObjects
     * @param   presentIndex            the present index of ParetoObject whose crowding distance is to be calculated
     * @param   objective               the objective function over which the value of i-th individual is to be calculated
     * @param   maxObjectiveValue       the maximum value for this objective function
     * @param   minObjectiveValue       the minimum value for this objective function
     * @return                          the crowding distance
     */
    private static double calculateCrowdingDistance(List<ParetoObject> singlePareto,
                                                    final int presentIndex,
                                                    final int objectiveNum,
                                                    final IObjectiveFunction objective,
                                                    final double maxObjectiveValue,
                                                    final double minObjectiveValue) {
        if(objectiveNum==0)
        return (
            singlePareto.get(presentIndex).getCrowdingDistance()
            + ((ObjFuction0.get(singlePareto.get(presentIndex + 1).getChromosome())
            - ObjFuction0.get(singlePareto.get(presentIndex - 1).getChromosome())) / (maxObjectiveValue - minObjectiveValue))
        );
        else return (
                singlePareto.get(presentIndex).getCrowdingDistance()
                + ((ObjFuction1.get(singlePareto.get(presentIndex + 1).getChromosome())
                - ObjFuction1.get(singlePareto.get(presentIndex - 1).getChromosome())) / (maxObjectiveValue - minObjectiveValue))
            ); 
    }
    
    /**
     * checks whether any of the dominated chromosome list of the given front is empty,
     * returns true if at least one set of dominated chromosomes is not non-empty.
     * 
     * @param   front   list of chromosomes whose dominated chromosomes are to be checked
     * @return          boolean logic whether the dominated chromosomes are empty
     */
    private static boolean isDominatedChromosomesEmpty(List<Chromosome> front) {
        return front.stream().anyMatch((chromosome) -> (!chromosome.getDominatedChromosomes().isEmpty()));
    }
    
    /**
     * checks if any of the dominated chromosomes of the input chromosome has a domination rank of 0,
     * returns true if at least one dominated chromosome contains domination rank 0.
     * 检查输入染色体的任何一个支配染色体的支配秩是否为0，如果至少有一个支配染色体的支配秩为0，则返回true
     * @param   chromosome  chromosome to check whether it contains any dominated chromosome with rank 0
     * @return  boolean logic whether dominated chromosomes contain rank 0.
     */
    private static boolean hasRecessiveRankGreaterThanZero(Chromosome chromosome) {
        
        if(chromosome.getDominatedChromosomes().isEmpty()) return false;
        
        return chromosome.getDominatedChromosomes().stream().noneMatch((recessive) -> (recessive.getDominationRank() == 0));
    }
    
    /**
     * the child and parent population is combined to create a larger population pool
     * 
     * @param   parent  parent population
     * @param   child   child population
     * @return          combined parent + child population
     */
    public static Population createCombinedPopulation(Population parent, Population child) {
        
        List<Chromosome> combinedPopulace = new ArrayList<>();
        Population combinedPopulation = new Population();

        combinedPopulace.addAll(parent.getPopulace());
        combinedPopulace.addAll(child.getPopulace());
        combinedPopulation.setPopulace(combinedPopulace);
        
        return combinedPopulation;
    }
    
    /**
     * this method decodes the genetic code that is represented as a string of binary values, converted into
     * decimal value.
     * 
     * @param   geneticCode     the genetic code as an array of Allele. Refer Allele.java for more information
     * @return                  the decimal value of the corresponding binary string.
     */
    public static double decodeGeneticCode(final Allele[] geneticCode) {//二进制换成十进制

        double value = 0;
        String binaryString = "";
        
        for(Allele bit : geneticCode) binaryString += bit.getGene() ? "1" : "0";
        for(int i = 0; i < binaryString.length(); i++) if(binaryString.charAt(i) == '1') value += Math.pow(2, binaryString.length() - 1 - i);
        
        return value;
    }
    
    /**
     * fitness is calculated using min-max normalization
     * 
     * @param   geneticCode     the genetic code whose fitness is to be calculated
     * @return                  the corresponding calculated fitness
     */
    public static double calculateFitness(Allele[] geneticCode) {
        return minMaxNormalization(decodeGeneticCode(geneticCode));
    }
    
    /**
     * an implementation of min-max normalization
     * 
     * @param   value   the value that is to be normalized
     * @return          the normalized value
     */
    private static double minMaxNormalization(final double value) {
        return (((value - Configuration.ACTUAL_MIN) / (Configuration.ACTUAL_MAX - Configuration.ACTUAL_MIN)) * (Configuration.NORMALIZED_MAX - Configuration.NORMALIZED_MIN)) + Configuration.NORMALIZED_MIN;
    }
    
    /**
     * used to generate a random integer value
     * 
     * @return a random integer value
     */
    public static int generateRandomInt() {
        return ThreadLocalRandom.current().nextInt();
    }
    
    /**
     * a short hand for System.out.println().
     * 
     * @param string    the string to print to console.
     */
    public static void p(String string) {
        System.out.println(string);
    }
    
    public static void clearObjFuction() {
    	ObjFuction0.clear();
        ObjFuction1.clear();
    }
    
    public static void writeFile(Population populace) throws IOException {
    	List<Chromosome> childPopulace=populace.getPopulace();
    	FileWriter fw = new FileWriter("D:\\推荐系统\\NSGA-Debug\\NSGA-II-master\\result.txt",true);	
    	fw.write("pareto最优集：" + "\r\n");
    	
    	for(Chromosome chromosome:childPopulace) {
    		double[] w = Synthesis.weights(chromosome.getGeneticCode());
    		
    		fw.write(chromosome.toString()+" ||\t "+w[0] + " "+w[1]+" || \t"+ObjFuction0.get(chromosome)+" "+ObjFuction1.get(chromosome)+"\r\n");
    	}
    	fw.write("\r\n");
    	fw.close();
    }
    
    public static double[] Topsis(Population populace) {
    	
    	List<Chromosome> childPopulace=populace.getPopulace();
    	double f0=0,f1=0;//作为标准化的分母
    	for(Chromosome chromosome:childPopulace) {
    		f0+=Math.pow(ObjFuction0.get(chromosome), 2);
    		f1+=Math.pow(ObjFuction1.get(chromosome), 2);
    	}
    	f0=Math.sqrt(f0);f1=Math.sqrt(f1);
    	
    	HashMap<Chromosome,Double[]> matrix1 = new HashMap<>();//加权标准化决策矩阵：染色体-目标函数加权标准化值
    	double v0max=Double.MIN_VALUE,v1max=Double.MIN_VALUE;//目标函数的负理想解(最大值)
    	double v0min=Double.MAX_VALUE,v1min=Double.MAX_VALUE;//目标函数的负理想解(最小值)
    	
    	for(Chromosome chromosome:childPopulace) {
    		double[] w=Synthesis.weights(chromosome.getGeneticCode());
    		Double[] array=new Double[2];
    		array[0]=(ObjFuction0.get(chromosome)/f0)*w[0];
    		if(v0max<array[0]) v0max=array[0];
    		if(v0min>array[0]) v0min=array[0];
    		array[1]=(ObjFuction1.get(chromosome)/f1)*w[1];
    		if(v1max<array[1]) v1max=array[1];
    		if(v1min>array[1]) v1min=array[1];
    		matrix1.put(chromosome, array);
    	}
    	
    	HashMap<Chromosome,Double[]> matrix2 = new HashMap<>();//染色体-到正理想解欧几里得距离distance[0]，负欧几里得距离distance[1]
    	for(Chromosome chromosome:childPopulace) {
    		Double[] distance=new Double[2]; 
    		distance[0]=Math.sqrt(Math.pow(matrix1.get(chromosome)[0]-v0max, 2)+Math.pow(matrix1.get(chromosome)[1]-v1max, 2));
    		distance[1]=Math.sqrt(Math.pow(matrix1.get(chromosome)[0]-v0min, 2)+Math.pow(matrix1.get(chromosome)[1]-v1min, 2));
    		matrix2.put(chromosome, distance);    		
    	}
    	
    	double cmax=Double.MIN_VALUE;//用来评出最优染色体
    	Chromosome priorChromosome=null;//存储最优值
    	for(Chromosome chromosome:childPopulace) {
    		double c = matrix2.get(chromosome)[1]/(matrix2.get(chromosome)[1]+matrix2.get(chromosome)[0]);
    		if(cmax<=c) {cmax=c;priorChromosome=chromosome;}
    	}
    	double[] priorWeight=Synthesis.weights(priorChromosome.getGeneticCode());
    	return  priorWeight;
	}
}

	





///*
// * This repository / codebase is Open Source and free for use and rewrite.
// */
//package io.onclave.nsga.ii.api;
////fastNonDominatedSort,calculateCrowdingDistance,calculateFitness
//import io.onclave.nsga.ii.Interface.IObjectiveFunction;
//import io.onclave.nsga.ii.configuration.Configuration;
//import io.onclave.nsga.ii.datastructure.Allele;
//import io.onclave.nsga.ii.datastructure.Chromosome;
//import io.onclave.nsga.ii.datastructure.ParetoObject;
//import io.onclave.nsga.ii.datastructure.Population;
//import net.librec.common.LibrecException;
//
//import java.io.FileWriter;
//import java.io.IOException;
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.HashMap;
//import java.util.List;
//import java.util.concurrent.ThreadLocalRandom;
//
////这是一个服务类，它完成了大部分底层工作，这些工作是由业务/控制器层的其他类抽象或封装的。
//
///**
// * This is the service class that does most of the under-the-hood work that is abstracted/encapsulated
// * by other classes at the business/controller layer.
// * 
// * @author  Debabrata Acharya <debabrata.acharya@icloud.com>
// * @version 1.1
// * @since   0.1
// */
//public class Service {
//    
//    /**
//     * this is an implementation of the fast non-dominated sorting algorithm as defined in the
//     * NSGA-II paper [DOI: 10.1109/4235.996017] Section III Part A.
//     * 
//     * @param   population  the population object that needs to undergo fast non-dominated sorting algorithm
//     * @return  a HashMap with an integer key that labels the ranks and a list of chromosomes as values that clot chromosomes of same rank
//     */
//	//三者中间量->245，为crowdingDistanceAssignment服务
//	static List<Chromosome>  chromToFuntionToCrowd = new ArrayList<>(); //存储染色体，与保存的目标函数值一一对应
//	static double[] crowdObjFuntion0;//存储所有染色体对应的函数1值ֵ
//	static double[] crowdObjFuntion1; //存储所有染色体对应的函数2值ֵ
//	
//    public static HashMap<Integer, List<Chromosome>> fastNonDominatedSort(Population population) {
//    	System.out.println("22222");
//    	
//        HashMap<Integer, List<Chromosome>> paretoFront = new HashMap<>();//F1，F2,F3......
//        List<Chromosome> singularFront = new ArrayList<>();//chromosome.getDominationRank() == 0
//        List<Chromosome> populace = population.getPopulace();
//        
//      //将所有解对应的目标函数值保存
//        System.out.println("求目标函数");
//        List<IObjectiveFunction> objectives = Configuration.getObjectives();
//        double[] chromObjFuntion0 = new double[2*Configuration.getPOPULATION_SIZE()+2];//存储所有染色体对应的函数1值ֵ
//        double[] chromObjFuntion1 = new double[2*Configuration.getPOPULATION_SIZE()+2];//存储所有染色体对应的函数2值ֵ
//        List<Chromosome>  chromToFuntion = new ArrayList<>(); 							//存储染色体，与保存的目标函数值一一对应
//        int count = 0;
//        for (Chromosome chromosome : populace) {
//        	//System.out.println(count);
//			chromObjFuntion0[count] = objectives.get(0).objectiveFunction(chromosome);
//			chromObjFuntion1[count] = objectives.get(1).objectiveFunction(chromosome);
//			chromToFuntion.add(chromosome);
//			count++;
//		}
//        System.out.println("求目标函数");
//        chromToFuntionToCrowd = chromToFuntion;
//        crowdObjFuntion0 = chromObjFuntion0;
//        crowdObjFuntion1 = chromObjFuntion1;
//        /**
//         * iterating over each chromosome of the population
//         *///System.out.println("1for-on");
//        for(Chromosome chromosome : populace) {
//        	//为每条染色体非配支配等级
//            /**
//             * an initial domination rank of 0 is set for each chromosome and a blank list is set for the number of
//             * chromosomes that the present chromosome dominates.
//             */
//            chromosome.setDominationRank(0);
//            chromosome.setDominatedChromosomes(new ArrayList<>());
//            /**
//             * for each chromosome, the program iterates over all the other remaining chromosomes to find which other
//             * chromosomes are dominated by this chromosome and vice versa.
//             *///System.out.println("2for-on");
//            System.out.println("33333");
//            for (Chromosome competitor : populace) 
//            	{if(!competitor.equals(chromosome)) {
//            	//求出每条染色体所支配的个体集
//                /**
//                 * if the present chromosome dominates the competitor, then:
//                 *      i:   check if the competitor already exists in the list of dominated chromosomes of the present chromosome.
//                 *     ii:   if the competitor does not exist within the list, then add it to the list of dominated chromosomes
//                 *           of the present chromosome.
//                 * else, if the competitor dominates the present chromosome, then increment the domination rank of the present
//                 * chromosome by one.
//                 *///System.out.println("domi-on");
//            		
//            	  	double[] chromosomeFunct = new double[2], competitorFunt = new double[2];   //存储查找到的目标函数值           	
//            	  	for(int i=0;i < 2*Configuration.getPOPULATION_SIZE();i++) {	//查找对应的chromosome目标函数值ֵ
//            	  		if(chromosome == chromToFuntion.get(i)) {
//            	  			chromosomeFunct[0] = chromObjFuntion0[i];
//            	  			chromosomeFunct[1] = chromObjFuntion1[i];break;
//            	  			}            	  		
//            	  	}
//            	  	
//            	  	for(int i=0;i < 2*Configuration.getPOPULATION_SIZE();i++) {	//查找对应的competitor目标函数值ֵ
//            	  		if(competitor == chromToFuntion.get(i)) {
//            	  			competitorFunt[0] = chromObjFuntion0[i];
//            	  			competitorFunt[1] = chromObjFuntion1[i];break;
//            	  			}            	  		
//            	  	}
//            	  	
//                if(dominates(chromosomeFunct, competitorFunt)) {//System.out.println("domi-over1");
//                    if(!chromosome.getDominatedChromosomes().contains(competitor)) chromosome.getDominatedChromosomes().add(competitor);
//                } 
//                else if(dominates(competitorFunt, chromosomeFunct)) {/*System.out.println("domi-over2");*/chromosome.setDominationRank(chromosome.getDominationRank() + 1);}
//                }}System.out.println("33333");
//            /**
//             * if the domination rank of the present chromosome is 0, it means that this chromosome is a non-dominated chromosome
//             * and hence it is added to the clot of chromosomes that are also non-dominated.
//             */
//            if(chromosome.getDominationRank() == 0) singularFront.add(chromosome);
//        }
//  //      System.out.println("1for-over");
//        /**
//         * the first clot of non-dominated chromosomes is added to the HashMap with rank label 1.
//         */
//        paretoFront.put(1, singularFront);
//        
//        int i = 1;
//        List<Chromosome> previousFront = paretoFront.get(i);
//        List<Chromosome> nextFront = new ArrayList<>();
//        
//        /**
//         * the current/previous ranked clot of chromosomes with rank i is iterated over to find the next clot of chromosomes
//         * with rank (i+1)
//         */
//        while(previousFront != null && !previousFront.isEmpty()) {
//            
//            /**
//             * iterating over each chromosome from the previous clot of chromosomes ranked i.
//             */
//            for(Chromosome chromosome : previousFront) {
//                
//                /**
//                 * iterating over each of the dominated chromosomes from the present chromosome of rank i.
//                 */
//                for(Chromosome recessive : chromosome.getDominatedChromosomes()) {
//                    
//                    /**
//                     * if the domination rank of the current recessive chromosome in consideration is not 0, then
//                     * decrement it's rank by 1.
//                     * if the domination rank of the current recessive chromosome in consideration is 0, then add
//                     * it to the next front [clot of chromosomes that belong to rank (i+1)].
//                     */
//                    if(recessive.getDominationRank() != 0) recessive.setDominationRank(recessive.getDominationRank() - 1);
//                    if(recessive.getDominationRank() == 0) if(!nextFront.contains(recessive)) nextFront.add(recessive);
//                }
//            }
//            
//            /**
//             * this code snippet ensures "rank jumps" to create all the possible rank lists from the parent
//             * population.
//             * new ranks are created only when there are recessive chromosomes with domination rank = 1 which are
//             * decremented to domination rank 0 and then added to the next front.
//             * but, due to the randomness of the algorithm, situation may occur such that even after decrementing all recessive
//             * chromosome domination ranks by 1, none have domination rank 0 and hence the next front remains empty.
//             * to ensure that all recessive chromosomes are added to some rank list, the program jumps domination ranks
//             * of each recessive chromosome by decrementing domination rank by 1 until at least one of them reaches a
//             * domination rank count of 0 and then that recessive chromosome is added to the next front.
//             * 
//             * if the next front is empty and the previous front has at least one dominated chromosome:
//             *      i:  find the minimum rank among all the recessive chromosomes available:
//             *              1:  iterate over all the chromosomes of the previous front
//             *              2:  while the chromosomes have no dominated chromosomes with rank 0:
//             *                      a:  iterate over all the recessive chromosomes of the current chromosome
//             *                      b:  if the minimum rank is greater than the dominated rank of the present recessive,
//             *                          mark this as the minimum rank recorded among all recessive chromosomes available.
//             *              3:  end while
//             *     ii:  iterate over all the chromosomes of the previous front
//             *              1: while the chromosomes have no dominated chromosomes with rank 0:
//             *                      a:  iterate over all the dominated chromosomes of the current chromosome
//             *                      b:  if the domination rank of the recessive chromosome is not 0, then decrement the
//             *                          domination count by value of minimum rank.
//             *                      c:  if the domination rank is 0, then add it to the next front.
//             *              2:  end while
//             */
//            if(nextFront.isEmpty() && !isDominatedChromosomesEmpty(previousFront)) {
//                
//                int minimumRank = -1;
//                
//                for(Chromosome chromosome : previousFront)
//                    while(hasRecessiveRankGreaterThanZero(chromosome))
//                        for(Chromosome recessive : chromosome.getDominatedChromosomes())
//                            if((minimumRank == -1) || minimumRank > recessive.getDominationRank()) minimumRank = recessive.getDominationRank();
//                
//                if(minimumRank != -1) for(Chromosome chromosome : previousFront)
//                    while(hasRecessiveRankGreaterThanZero(chromosome)) for(Chromosome recessive : chromosome.getDominatedChromosomes()) {
//                            if(recessive.getDominationRank() != 0) recessive.setDominationRank(recessive.getDominationRank() - minimumRank);
//                            if(recessive.getDominationRank() == 0) if(!nextFront.contains(recessive)) nextFront.add(recessive);
//                    }
//            }
//            
//            /**
//             * if the next front calculated is not empty, then it is added to the ranked HashMap data-structure
//             * with the rank (i+1), else all chromosomes are sorted into some rank or the other and the program
//             * breaks out of the loop.
//             */
//            if(!nextFront.isEmpty()) paretoFront.put(++i, nextFront); else break;
//            
//            /**
//             * the next front (i) calculated is marked as the previous front for the next iteration (i+1) and
//             * an empty next front is created.
//             */
//            previousFront = nextFront;
//            nextFront = new ArrayList<>();
//        }
//        
//        return paretoFront;
//    }
//    
//    /**
//     * this is the implementation of the crowding distance assignment algorithm as defined in the
//     * NSGA-II paper [DOI: 10.1109/4235.996017] Section III Part B.
//     * this ensures diversity preservation.
//     * 
//     * @param   singularFront   a list of chromosomes whose crowding distances are to be calculated
//     * @return                  a list of ParetoObjects with assigned crowding distances. [Refer ParetoObject.java for more information]
//     * @throws LibrecException 
//     */
//    public static List<ParetoObject> crowdingDistanceAssignment(List<Chromosome> singularFront) {
//        
//        int i = 0;
//        int end = singularFront.size() - 1;
//        Double maxObjectiveValue;
//        Double minObjectiveValue;
//        List<IObjectiveFunction> objectives = Configuration.getObjectives();
//        List<ParetoObject> singlePareto = new ArrayList<>();
//        
//        /**
//         * for each chromosome in the input list, a new ParetoObject with an initial crowding distance of 0
//         * is created and added to the list of ParetoObjects that are to be returned.
//         */
//        for(Chromosome chromosome : singularFront) singlePareto.add(i++, new ParetoObject(chromosome, 0f));
//        //ӵ���ȳ�ʼ��Ϊ0
//        /**
//         * iterating over each of the objective functions set [refer Configuration.java for more information],
//         * the ParetoObject list is sorted according to the objective functions and the first and last ParetoObjects
//         * are set a crowding distance of infinity.
//         */ int objectiveNum = 0;
//        for(IObjectiveFunction objective : objectives) {
//           
//           double[] paretoObjFuntionNum;
//           if(objectiveNum == 0) paretoObjFuntionNum = crowdObjFuntion0;	//存储所有染色体对应的函数1值ֵ
//           else paretoObjFuntionNum = crowdObjFuntion1;						//存储所有染色体对应的函数2值ֵ
//           objectiveNum++;
//           
//            maxObjectiveValue = null;
//            minObjectiveValue = null;
//            singlePareto = sort(singlePareto, objective);
//            
//            singlePareto.get(0).setCrowdingDistance(Double.MAX_VALUE);
//            singlePareto.get(end).setCrowdingDistance(Double.MAX_VALUE);
//            
//            /**
//             * the max and min objective values are calculated according to the present objective function
//             *///System.out.println("cr-o");
//            double[] singleParetoFunct = new double[singlePareto.size()];
//            for(ParetoObject paretoObject : singlePareto) {
//            	          	
//            	double paretoObjectFunct = 0;   //存储查找到的目标函数值          	
//        	  	for(int k=0;k < 2*Configuration.getPOPULATION_SIZE();k++) {	//查找对应的paretoObject目标函数值ֵ
//        	  		if(paretoObject.getChromosome() == chromToFuntionToCrowd.get(k)) {
//        	  			paretoObjectFunct = paretoObjFuntionNum[k]; break;
//        	  			}            	  		
//        	  	}
//        	  //找到paretoObject与目标值一一对应
//        	  	for(int k=0,j = 0;k < 2*Configuration.getPOPULATION_SIZE();k++) {	//查找对应的chromosome目标函数值ֵ
//        	  		if(paretoObject.getChromosome() == chromToFuntionToCrowd.get(k)) {
//        	  			singleParetoFunct[j] = paretoObjFuntionNum[k]; j++;
//        	  			}            	  		
//        	  	}
//        	  	
//                if((maxObjectiveValue == null) || (maxObjectiveValue < paretoObjectFunct)) maxObjectiveValue = paretoObjectFunct;
//                if((minObjectiveValue == null) || (minObjectiveValue > paretoObjectFunct)) minObjectiveValue = paretoObjectFunct;
//            }
//           // System.out.println("cr-end");
//            /**
//             * the crowding distance of all ParetoObjects are calculated and assigned except the first and last ParetoObjects
//             * that have infinite crowding distance
//             */
//            for(i = 2; i < end; i++) singlePareto.get(i).setCrowdingDistance(calculateCrowdingDistance(singlePareto,i,objective,maxObjectiveValue,minObjectiveValue,singleParetoFunct));
//        }
//        //System.out.println("crowdingend");
//        return singlePareto; 
//   }
//    
//    /**
//     * this method sorts a list of ParetoObjects based on the Crowd-Comparison Operator using the domination
//     * rank and crowding distance as discussed in the NSGA-II paper [DOI: 10.1109/4235.996017] Section III Part B.
//     * 
//     * @param   singleFront     a list of ParetoObjects that are to be sorted according to their crowding distance
//     * @return                  a list of sorted ParetoObjects
//     */
//    public static List<ParetoObject> crowdComparisonSort(List<ParetoObject> singleFront) {
//        
//        int index = -1;
//        List<ParetoObject> sortedFront = new ArrayList<>();
//        ParetoObject presentParetoObject;
//        ParetoObject competitor;
//        
//        /**
//         * all the ParetoObjects are, at first, marked as false for crowding distance sorted.
//         */
//        singleFront.stream().forEach((paretoObject) -> { paretoObject.setCrowdingDistanceSorted(false); });
//        
//        /**
//         * iterating over each ParetoObject in the singular front input:
//         *  i:  the i-th ParetoObject is marked as presentParetoObject
//         * ii:  if the presentParetoObject is not already sorted by crowding distance:
//         *          1:  iterate over the rest of the ParetoObjects in the input list as competitors that are
//         *              not already sorted using crowding distance
//         *          2:  compare the i-th and the j-th chromosome using the crowd comparison operator:
//         *                  a: for different ranks, choose the one with the lower (better) rank.
//         *                  b: for same rank, choose the one which has lower crowding distance.
//         *          3:  if competitor dominates the i-th chromosome, then mark competitor as presentParetoObject
//         *          4:  continue until i-th chromosome is compared to all competitors.
//         *          5:  mark the presentParetoObject as already sorted by crowding distance
//         *          6:  add presentParetoObject into list of sorted front with an incremented index
//         */
//        for(int i = 0; i < singleFront.size(); i++) {
//            
//            presentParetoObject = singleFront.get(i);
//            
//            if(!presentParetoObject.isCrowdingDistanceSorted()) {
//                
//                for(int j = 0; j < singleFront.size(); j++) {
//
//                    competitor = singleFront.get(j);
//                    
//                    if(!competitor.isCrowdingDistanceSorted()) {
//                        
//                        double dominationRank = presentParetoObject.getChromosome().getDominationRank();
//                        double competingDominationRank = competitor.getChromosome().getDominationRank();
//                        double crowdingDistance = presentParetoObject.getCrowdingDistance();
//                        double competingCrowdingDistance = competitor.getCrowdingDistance();
//
//                        if(i != j) if((dominationRank > competingDominationRank) || ((dominationRank == competingDominationRank) && (crowdingDistance < competingCrowdingDistance))) presentParetoObject = competitor;
//                    }
//                }
//                
//                presentParetoObject.setCrowdingDistanceSorted(true);
//                sortedFront.add(++index, presentParetoObject);
//            }
//        }
//        
//        return sortedFront;
//    }
//    
//    /**
//     * this method is not implemented, as it is not absolutely necessary for this algorithm to work.
//     * is kept if implementation is needed in future.
//     * returns the same unsorted parent population as of now.
//     * 
//     * @param   population  the population that is to be sorted
//     * @return              a sorted population
//     */
//    public static Population nonDominatedPopulationSort(Population population) {
//        
//        //--TO-DO--
//        
//        return population;
//    }
//    
//    /**
//     * this method checks whether competitor1 dominates competitor2.
//     * requires that none of the values of the objective functions using competitor1 is smaller
//     * than the values of the objective functions using competitor2.
//     * at least one of the values of the objective functions using competitor1 is greater than
//     * the corresponding value of the objective functions using competitor2.
//     * 
//     * @param   competitor1     the chromosome that may dominate
//     * @param   competitor2     the chromosome that may be dominated
//     * @return                  boolean logic whether competitor1 dominates competitor2.
//     */
//    public static boolean dominates(final double[] competitorFunt1, final double[] competitorFunt2) {
// //   	System.out.println("iner-on");
//    	//long time1 = System.currentTimeMillis();
//        /**
//         * getting the list of configured objectives from Configuration.java
//         */
//    	if(competitorFunt1[0] > competitorFunt2[0] && competitorFunt1[1] > competitorFunt2[1]) return true;
//    	else return false;
//    	
//      /*  List<IObjectiveFunction> objectives = Configuration.getObjectives();
//        double F1 = objectives.get(0).objectiveFunction(competitor1);
//        double G1 = objectives.get(1).objectiveFunction(competitor1);
//        double F2 = objectives.get(0).objectiveFunction(competitor2);
//        double G2 = objectives.get(1).objectiveFunction(competitor2);*/
//        
//        //long time2 = System.currentTimeMillis();
//        
//       // if(F1 > F2 && G1 > G2) 
//       	/*{System.out.println("比较一次时间："+(time2-time1));
//        	System.out.println("if");System.out.println(F1+"\t"+F2+"\n"+G1+"\t"+G2);System.out.println(System.currentTimeMillis());
//        	return true;}*/
//       /* else 
//        	{System.out.println("比较一次时间："+(time2-time1));
//        	System.out.println("else");System.out.println(F1+"\t"+F2+"\n"+G1+"\t"+G2);System.out.println(System.currentTimeMillis());
//        	return false;}*��
//        /**
//         * checks the negation of the predicate [none of the values of objective functions using competitor1
//         * is less than values of objective functions using competitor2] meaning that at least one of the values
//         * of the objective functions using competitor1 is less than the values of the objective functions using
//         * competitor2, hence returning false as competitor1 does not dominate competitor2
//         *///false = 至少有一个目标函数(competitor1) < 目标函数(competitor2)
//       /*if (!objectives.stream().noneMatch((objective) -> (objective.objectiveFunction(competitor1) < objective.objectiveFunction(competitor2)))) 
//    	   {System.out.println("iner1-over");
//    	   System.out.println(objectives.get(0).objectiveFunction(competitor1)+"\t"+objectives.get(0).objectiveFunction(competitor2));
//    	   System.out.println(objectives.get(1).objectiveFunction(competitor1)+"\t"+objectives.get(1).objectiveFunction(competitor2));
//    	   System.out.println(System.currentTimeMillis());return false;}  
//       	
//        /** 
//         * returns the value of the predicate [at least one of the values of the objective functions using
//         * competitor1 is greater than the corresponding value of the objective function using competitor2]
//         *///true = 所有的目标函数(competitor1) > 目标函数(competitor2)
//       /*System.out.println("iner2-over");
//       System.out.println(objectives.get(0).objectiveFunction(competitor1)+"\t"+objectives.get(0).objectiveFunction(competitor2));
//       System.out.println(objectives.get(1).objectiveFunction(competitor1)+"\t"+objectives.get(1).objectiveFunction(competitor2));
//       boolean b = objectives.stream().anyMatch((objective) -> (objective.objectiveFunction(competitor1) > objective.objectiveFunction(competitor2)));
//       //System.out.println(b);
//       System.out.println(System.currentTimeMillis());
//       return objectives.stream().anyMatch((objective) -> (objective.objectiveFunction(competitor1) > objective.objectiveFunction(competitor2)));*/       	
//    }
//    
//    /**
//     * the list is first converted to an array data-structure and then a randomized quick sort
//     * algorithm is followed.
//     * the resulting sorted array is again converted to a List data-structure before returning.
//     * 
//     * @param   singlePareto    the list of ParetoObjects that are to be sorted.
//     * @param   objective       the objective function using which the ParetoObjects are sorted.
//     * @return                  sorted list of ParetoObjects.
//     * @throws LibrecException 
//     */
//    private static List<ParetoObject> sort(List<ParetoObject> singlePareto, IObjectiveFunction objective) {
//        
//        ParetoObject[] paretoArray = new ParetoObject[singlePareto.size()];
//        singlePareto.toArray(paretoArray);
//        
//        randomizedQuickSort(paretoArray, 0, paretoArray.length - 1, objective);
//        
//        return (new ArrayList<>(Arrays.asList(paretoArray)));
//    }
//    
//    /**
//     * refer [https://jordanspencerwu.github.io/randomized-quick-sort/] for more details on randomized
//     * quick sort algorithm.
//     * 
//     * @param   paretoArray     the array to be sorted
//     * @param   head            the pointer/position of the head element
//     * @param   tail            the pointer/position of the tail element
//     * @param   objective       the objective function depending on which the sort is to take place
//     * @return                  the pivot index.
//     * @throws LibrecException 
//     */
//    private static int partition(ParetoObject[] paretoArray, int head, int tail, IObjectiveFunction objective) {
//        
//        ParetoObject pivot = paretoArray[tail];
//        int i = head - 1;
//        
//        for(int j = head; j <= (tail - 1); j++) {
//            
//            if(objective.objectiveFunction(paretoArray[j]) <= objective.objectiveFunction(pivot)) {
//                
//                i++;
//                ParetoObject temporary = paretoArray[i];
//                paretoArray[i] = paretoArray[j];
//                paretoArray[j] = temporary;
//            }
//        }
//        
//        ParetoObject temporary = paretoArray[i + 1];
//        paretoArray[i + 1] = paretoArray[tail];
//        paretoArray[tail] = temporary;
//        
//        return (i + 1);
//    }
//    
//    /**
//     * refer [https://jordanspencerwu.github.io/randomized-quick-sort/] for more details on randomized
//     * quick sort algorithm.
//     * 
//     * @param   paretoArray     the array to be sorted
//     * @param   head            the pointer/position of the head element
//     * @param   tail            the pointer/position of the tail element
//     * @param   objective       the objective function depending on which the sort is to take place
//     * @return                  the random partition position index.
//     * @throws LibrecException 
//     */
//    private static int randomizedPartition(ParetoObject[] paretoArray, int head, int tail, IObjectiveFunction objective) {
//        
//        int random = ThreadLocalRandom.current().nextInt(head, tail + 1);
//        
//        ParetoObject temporary = paretoArray[head];
//        paretoArray[head] = paretoArray[random];
//        paretoArray[random] = temporary;
//        
//        return partition(paretoArray, head, tail, objective);
//    }
//    
//    /**
//     * refer [https://jordanspencerwu.github.io/randomized-quick-sort/] for more details on randomized
//     * quick sort algorithm.
//     * 
//     * @param   paretoArray     the array to be sorted
//     * @param   head            the pointer/position of the head element
//     * @param   tail            the pointer/position of the tail element
//     * @param   objective       the objective function depending on which the sort is to take place
//     * @throws LibrecException 
//     */
//    private static void randomizedQuickSort(ParetoObject[] paretoArray, int head, int tail, IObjectiveFunction objective) {
//        
//        if(tail < head) {
//            
//            int pivot = randomizedPartition(paretoArray, head, tail, objective);
//            
//            randomizedQuickSort(paretoArray, head, pivot - 1, objective);
//            randomizedQuickSort(paretoArray, pivot + 1, tail, objective);
//        }
//    }
//    
//    /**
//     * implementation of crowding distance calculation as defined in NSGA-II paper
//     * [DOI: 10.1109/4235.996017] Section III Part B.
//     * 
//     * I[i]distance = I[i]distance + (I[i+1].m - I[i-1].m)/(f-max - f-min)
//     * 
//     * I[i]distance = crowding distance of the i-th individual
//     * I[i+1].m = m-th objective function value of the (i+1)-th individual
//     * I[i-1].m = m-th objective function value of the (i-1)-th individual
//     * f-max, f-min = maximum and minimum values of the m-th objective function
//     * 
//     * @param   singlePareto            the list of ParetoObjects
//     * @param   presentIndex            the present index of ParetoObject whose crowding distance is to be calculated
//     * @param   objective               the objective function over which the value of i-th individual is to be calculated
//     * @param   maxObjectiveValue       the maximum value for this objective function
//     * @param   minObjectiveValue       the minimum value for this objective function
//     * @return                          the crowding distance
//     * @throws LibrecException 
//     */
//    private static double calculateCrowdingDistance(List<ParetoObject> singlePareto,
//                                                    final int presentIndex,
//                                                    final IObjectiveFunction objective,
//                                                    final double maxObjectiveValue,
//                                                    final double minObjectiveValue,
//                                                    double[] singleParetoFunct) {
//        
//        return (
//            singlePareto.get(presentIndex).getCrowdingDistance()
//            + ((singleParetoFunct[presentIndex + 1]
//            - singleParetoFunct[presentIndex - 1]) / (maxObjectiveValue - minObjectiveValue))
//        );
//    }
//    
//    /**
//     * checks whether any of the dominated chromosome list of the given front is empty,
//     * returns true if at least one set of dominated chromosomes is not non-empty.
//     * 
//     * @param   front   list of chromosomes whose dominated chromosomes are to be checked
//     * @return          boolean logic whether the dominated chromosomes are empty
//     */
//    private static boolean isDominatedChromosomesEmpty(List<Chromosome> front) {
//        return front.stream().anyMatch((chromosome) -> (!chromosome.getDominatedChromosomes().isEmpty()));
//    }
//    
//    /**
//     * checks if any of the dominated chromosomes of the input chromosome has a domination rank of 0,
//     * returns true if at least one dominated chromosome contains domination rank 0.
//     * 
//     * @param   chromosome  chromosome to check whether it contains any dominated chromosome with rank 0
//     * @return  boolean logic whether dominated chromosomes contain rank 0.
//     */
//    private static boolean hasRecessiveRankGreaterThanZero(Chromosome chromosome) {
//        
//        if(chromosome.getDominatedChromosomes().isEmpty()) return false;
//        
//        return chromosome.getDominatedChromosomes().stream().noneMatch((recessive) -> (recessive.getDominationRank() == 0));
//    }
//    
//    /**
//     * the child and parent population is combined to create a larger population pool
//     * 
//     * @param   parent  parent population
//     * @param   child   child population
//     * @return          combined parent + child population
//     */
//    public static Population createCombinedPopulation(Population parent, Population child) {
//        List<Chromosome> combinedPopulace = new ArrayList<>();
//        Population combinedPopulation = new Population();
//
//        combinedPopulace.addAll(parent.getPopulace());
//        combinedPopulace.addAll(child.getPopulace());
//        combinedPopulation.setPopulace(combinedPopulace);
//        return combinedPopulation;
//    }
//    
//    /**
//     * this method decodes the genetic code that is represented as a string of binary values, converted into
//     * decimal value.
//     * 
//     * @param   geneticCode     the genetic code as an array of Allele. Refer Allele.java for more information
//     * @return                  the decimal value of the corresponding binary string.
//     */
//    
//    public static double decodeGeneticCode(final Allele[] geneticCode) {//geneticCodeNumbers=10bit
//        double value = 0;
//        String binaryString = "";
//        
//        for(Allele bit : geneticCode) binaryString += bit.getGene() ? "1" : "0";
//        for(int i = 0; i < binaryString.length(); i++) if(binaryString.charAt(i) == '1') value += Math.pow(2, binaryString.length() - 1 - i);
//      //二进制基因转换成是十进制
//        value=value/Math.pow(2, binaryString.length());//产生value=[0,1]->weightArray
//        value=(double)Math.round(value*100)/100;//double后值保留两位->weightArray
//        return value;
//    }
//    
//  //产生权值得方法，调用时记得改方法decodeGeneticCode里的value值ֵ
//    /*static double[] Weights=new double[3];
//    public static void weightArray(){
//    	do{
//    		Weights[0]=decodeGeneticCode(Synthesis.synthesizeGeneticCode(Configuration.getCHROMOSOME_LENGTH()));
//    		Weights[1]=decodeGeneticCode(Synthesis.synthesizeGeneticCode(Configuration.getCHROMOSOME_LENGTH()));
//    		Weights[2]=decodeGeneticCode(Synthesis.synthesizeGeneticCode(Configuration.getCHROMOSOME_LENGTH()));
//    	}while((Weights[0]+Weights[1]+Weights[2])!=1);
//    }*/
//    
//    
//    /**
//     * fitness is calculated using min-max normalization
//     * 
//     * @param   geneticCode     the genetic code whose fitness is to be calculated
//     * @return                  the corresponding calculated fitness
//     */
//    public static double calculateFitness(Allele[] geneticCode) {
//        return minMaxNormalization(decodeGeneticCode(geneticCode));
//    }
//    
//    /**
//     * an implementation of min-max normalization
//     * 
//     * @param   value   the value that is to be normalized
//     * @return          the normalized value
//     */
//    private static double minMaxNormalization(final double value) {//value=(10 1010 1010)D
//        return (((value - Configuration.ACTUAL_MIN) / (Configuration.ACTUAL_MAX - Configuration.ACTUAL_MIN)) * (Configuration.NORMALIZED_MAX - Configuration.NORMALIZED_MIN)) + Configuration.NORMALIZED_MIN;
//    }//fitness function
//    
//    /**
//     * used to generate a random integer value
//     * 
//     * @return a random integer value
//     */
//    public static int generateRandomInt() {
//        return ThreadLocalRandom.current().nextInt();
//    }
//    
//    /**
//     * a short hand for System.out.println().
//     * 
//     * @param string    the string to print to console.
//     */
//    public static void p(String string) {
//        System.out.println(string);
//    }
//    
//    public static void writeFile(List<Chromosome> childPopulace,int k) throws IOException {
//    	FileWriter fw = new FileWriter("D:\\推荐系统\\NSGA-Debug\\NSGA-II-master\\result.txt",true);	
//    	fw.write("GENERATION : " + k + "\r\n");
//    	
//    	for(int j = 0;j < Configuration.getPOPULATION_SIZE();j++) {
//    		for(int i = 0;i < chromToFuntionToCrowd.size();i++) {
//    			if(childPopulace.get(j) == chromToFuntionToCrowd.get(i)) {
//    				double[] w = Synthesis.weights(childPopulace.get(j).getGeneticCode());   				
//					fw.write(w[0] + " "+w[1]+"\t\t"+crowdObjFuntion0[i]+" "+crowdObjFuntion1[i]+"\r\n");
//    			}
//    		}
//    	}
//    	fw.write("\r\n");
//    	fw.close();
//    }
//    
//    //TestWeights
//    public static void main(String[] args){//测试时记得改方法decodeGeneticCode里的value值ֵ
//
//    }
//}
