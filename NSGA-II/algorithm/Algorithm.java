/*
 * This repository / codebase is Open Source and free for use and rewrite.
 */
package io.onclave.nsga.ii.algorithm;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
//import org.jfree.ui.RefineryUtilities;

import org.jfree.ui.RefineryUtilities;

import Demo.Calculator;
import io.onclave.nsga.ii.api.GraphPlot;
import io.onclave.nsga.ii.api.Reporter;
import io.onclave.nsga.ii.api.Service;
import io.onclave.nsga.ii.api.Synthesis;
import io.onclave.nsga.ii.configuration.Configuration;
import io.onclave.nsga.ii.datastructure.Chromosome;
import io.onclave.nsga.ii.datastructure.ParetoObject;
import io.onclave.nsga.ii.datastructure.Population;
import net.librec.common.LibrecException;
import net.librec.recommender.AbstractRecommender;
import Demo.RecAlgorithm; 

/**
 * This is the starting point of the main NSGA-II algorithm.
 * Run this class to get the desired output.
 * 
 * @author  Debabrata Acharya <debabrata.acharya@icloud.com>
 * @version 1.1
 * @since   0.1
 */
public class Algorithm {
    
    /**
     * This method first prepares the multi-objectives that it is going to work with. In this case,
     * it works with 2 objectives. At generation 0, a random initial population is generated and then
     * sorted using non-dominated population sorting. Using this initial parent population, a child
     * population is generated. As the initial parent and child population are created, new generations
     * are simulated and at each generation, the following actions are carried out:
     *      1:  the present parent and child are combined to create a new population containing all
     *          chromosomes from both parent and child. This ensures elitism.
     *      2:  this new combined population is then sorted using the fast non-dominated sorting algorithm
     *          to get a list of chromosomes that are grouped according to their rank. The higher the rank,
     *          more desirable they are to be carried forward into the next generation.
     *      3:  an iteration is carried over all the ranks as follows:
     *              i:  the list of chromosomes from the current iterated rank is taken into account.
     *             ii:  the amount of free remaining space in the new child population is calculated.
     *            iii:  a crowd comparison sort is done after assigning crowding distance to the chromosomes.
     *             iv:  if the number of chromosomes in this rank is less than or equal to the amount of free
     *                  space available in the new child population, then the whole chromosome cluster is added
     *                  to the new population, else only the available number of chromosomes is added to the
     *                  child population according to their crowding distance. This is done for diversity
     *                  preservation.
     *              v:  end.
     *      4:  the new synthesized populace is added to the new child population.
     *      5:  if this is the last generation, then the present child is shown as the Pareto Front, otherwise,
     *          the present child is labeled as the new parent for the next generation, and a new child
     *          population is generated from this newly labeled parent population. This combination now becomes
     *          the new parent/child for the next generation.
     *      6:  all the child from all the generations are added to the Graph Rendering engine to show all the
     *          child data as fronts for that generation.
     *      7:  end.
     * the plotted graphs are viewed.
     * 
     * @param   args    pass command line arguments. Not required to run this code.
     * @throws LibrecException 
     * @throws IOException 
     * @throws NoSuchMethodException 
     * @throws InvocationTargetException 
     * @throws IllegalAccessException 
     * @throws ClassNotFoundException 
     * @see             Plotted graphs of all fronts as well as the Pareto Front as output. 
     */
	
    public static void main(String[] args) throws LibrecException, IOException, ClassNotFoundException, IllegalAccessException, InvocationTargetException, NoSuchMethodException {
    	
    	long time1 = System.currentTimeMillis();
    	
    	
    	RecAlgorithm.RecBuildModel();
//		Calculator.CalculatorRec(); 		
    	
        /* prepares the objectives [See Configuration.java file for more information.] */
       Configuration.buildObjectives();
       GraphPlot multiPlotGraph = new GraphPlot();
       
        /**
         * a new random population is synthesized and sorted using non-dominated population sort to get
         * a sorted list of parent chromosomes at generation 0.
         * child population generated from parent population.
         */
        Population parent = Service.nonDominatedPopulationSort(Synthesis.syntesizePopulation()); //设置了所有属性
        Population child = Synthesis.synthesizeChild(parent);//选择、交叉、变异生成下一代
       /* for(int k=0;k<150;k++)
		{double[] weights = Synthesis.weights(child.getPopulace().get(k).getGeneticCode());
		System.out.println(weights[0]+" "+weights[1]+"\t=\t"+(weights[0]+weights[1])+"\n");
		}*/
    	
        /**
         * a loop is run that iterates as new generations are created (new child population is created from previous parent
         * population.
         * the number of generations to be simulated are defined in the Configuration.java file.
         */
        for(int i = 2; i <= Configuration.getGENERATIONS(); i++) {
            
            System.out.println("GENERATION : " + i);
            
            /**
             * a combined population of both latest parent and child is created to ensure elitism.
             * the combined population created is then sorted using fast non-dominated sorting algorithm,
             * to create rank wise divisions [chromosomes with rank 1 (non-dominated),
             * chromosomes with rank 2 (dominated by 1 chromosome), etc.]
             * this information is stored in a HashMap data-structure that maps one integer value
             * to one list of chromosomes. The integer refers to the rank number while the list refers
             * to the chromosomes that belong to that rank.
             *///System.out.println("fast-on");
            long T1 = System.currentTimeMillis();
            HashMap<Integer, List<Chromosome>> rankedFronts = Service.fastNonDominatedSort(Service.createCombinedPopulation(parent, child));
            long T2 = System.currentTimeMillis();System.out.println((T2-T1)/1000+"s");
            //染色体集合F1,F2,F3......且合成种群大小2N（已非支配排序）			//调用了目标函数
         //System.out.println("fast-over");
            Population nextChildPopulation = new Population();//最优的，含有childPopulace
            List<Chromosome> childPopulace = new ArrayList<>();//从rankedFronts（F1,F2,F3,...）中选出前N个子代   
            /**
             * an iteration is carried over the HashMap to go through each rank of chromosomes, and the
             * most desired chromosomes (higher ranks) are included into the child population of the
             * next generation.
             */
            for(int j = 1; j <= rankedFronts.size(); j++) {
                //从(最优非支配集合)F1,F2,F3...中选出最优的前N个染色体给childPopulace
                /**
                 * during iteration, the current ranked list of chromosomes is chosen and the amount of
                 * free space (to accommodate chromosomes) of the current child population is calculated
                 * to check whether chromosomes from this rank can be fit into the new child population.
                 */
            	
                List<Chromosome> singularFront = rankedFronts.get(j);
                int usableSpace = Configuration.getPOPULATION_SIZE() - childPopulace.size();
             
                /**
                 * if the new list of chromosomes is not null and if the child population has free usable space,
                 * then an attempt to include some or all of the chromosomes is made otherwise, there is no more
                 * space in the child population and hence no more rank/chromosome checks are done and the program
                 * breaks out from the inner for-loop.
                 */
                //System.out.println("usableSpace-on");
                
                if(singularFront != null && !singularFront.isEmpty() && usableSpace > 0) {
                
                    /**
                     * if the amount of usable space is more than or equal to the number of chromosomes in the clot,
                     * the whole clot of chromosomes is added to the child population/populace, otherwise, only the
                     * number of chromosomes that can be fit within the usable space is chosen according to the
                     * crowding distance of the chromosomes.
                     */
                    if(usableSpace >= singularFront.size()) childPopulace.addAll(singularFront);
                    else {
                        //当解集Fi中大小>N时的操作，拥挤距离排序，从解集Fi中选出解个体，使得子代种群childPopulace大小为N
                        /**
                         * a crowd comparison sort is carried over the present clot of chromosomes after assigning them a
                         * crowding distance (to preserve diversity) and hence a list of ParetoObjects are prepared.
                         * [refer ParetoObject.java for more information]
                         */
                    	//System.out.println("crowd-on");singularFront=解集Fi
                        List<ParetoObject> latestFront = Service.crowdComparisonSort(Service.crowdingDistanceAssignment(singularFront));
                        //System.out.println("crowd-over");									//调用目标函数
                        for(int k = 0; k < usableSpace; k++) 
//                        	if(!childPopulace.contains(latestFront.get(k).getChromosome()))
                        		childPopulace.add(latestFront.get(k).getChromosome());
                    }
                } else break;
                
            }
            //System.out.println("usableSpace-over");
            /**
             * the new populace is added to the new child population
             */
            nextChildPopulation.setPopulace(childPopulace);//精英策略：本次迭代由结合种群最终产生的能够繁衍的子代种群，不是下一代，是用来生成下一代种群
														//下面的种群child即为本代种群
//            System.out.println(childPopulace.size());
            //Service.writeFile(childPopulace,i);
            
            /**
             * if this iteration is not the last generation, the new child created is made the parent for the next
             * generation, and a new child is synthesized from this new parent for the next generation.
             * this is the new parent and child for the next generation.
             * if this is the last generation, no new parent/child combination is created, instead the Pareto Front
             * is plotted and rendered as the latest created child is the actual Pareto Front.
             */
            if(i < Configuration.getGENERATIONS()) {
                parent = child;           
                child = Synthesis.synthesizeChild(nextChildPopulation);  
            	} 
            else {
            	Reporter.render2DGraph(child);
            	Service.writeFile(child);//child即为本代种群
            	
            	//下面是进行TOPSIS算法
            	double[] priorWeight=Service.Topsis(child);
                System.out.println(priorWeight[0]+" "+priorWeight[1]);
                AbstractRecommender recommender=RecAlgorithm.CalculatorItemList(priorWeight);
                System.out.println("Percision:" + RecAlgorithm.getPrecision(recommender)*18);
                System.out.println("Recall:" + RecAlgorithm.getRecall(recommender));
//                AbstractRecommender recommender = Calculator.CalculatorItemList(priorWeight);
//                System.out.println("Percision:" + Calculator.getPercision(recommender));
//                System.out.println("Recall:" + Calculator.getReciprocal(recommender));
            	}
//本次迭代过程结束         
            Service.clearObjFuction();//迭代完清空本次集合，节约存储

            /**
             * this adds the child of each generation to the plotting to render the front of all the generations.
             */
            multiPlotGraph.prepareMultipleDataset(child, i, "generation " + i);
            
            //创建一个矩阵Z[][]存储每一个目标函数值ֵ
//            if( i == Configuration.getGENERATIONS()) {
//            	for(int k = 0;k < Configuration.getPOPULATION_SIZE();k++)
//            		for(int j = 0;j < Configuration.getObjectives().size();j++) {
//            			Z[k][j] = Configuration.getObjectives().get(j).objectiveFunction(child.getPopulace().get(k));
//            		}
//
//                //输出最优权值ֵ
//                int[] c = Synthesis.Topsis(Z);
//                double[] weights = Synthesis.weights(child.getPopulace().get(c[0]).getGeneticCode());
//                System.out.println(weights[0]+" "+weights[1]+"\t=\t"+(weights[0]+weights[1])+"\n");
//                
//                AbstractRecommender recommender = Calculator.CalculatorItemList(weights);
//                System.out.println("Percision:" + Calculator.getPercision(recommender));
//                System.out.println("Reciprocal:" + Calculator.getReciprocal(recommender));
//            }
            
        }
        
        System.out.println("\n\n----CHECK PARETO FRONT OUTPUT----\n\n");
        
        long time2 = System.currentTimeMillis();
        
        System.out.println((time2-time1)/1000 + "s");
        
        /**
         * the plotted and rendered chart/graph is viewed to the user.
         */
       multiPlotGraph.configureMultiplePlotter(Configuration.getXaxisTitle(), Configuration.getYaxisTitle(), "All Pareto");
       multiPlotGraph.pack();
       RefineryUtilities.centerFrameOnScreen(multiPlotGraph);
       multiPlotGraph.setVisible(true);
    }
}
