package io.onclave.nsga.ii.objectivefunction;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;

import Demo.Calculator;
import Demo.RecAlgorithm;
import io.onclave.nsga.ii.Interface.IObjectiveFunction;
import io.onclave.nsga.ii.api.Synthesis;
import io.onclave.nsga.ii.datastructure.Allele;
import io.onclave.nsga.ii.datastructure.Chromosome;
import io.onclave.nsga.ii.datastructure.ParetoObject;
import net.librec.common.LibrecException;
import net.librec.recommender.*;
import net.librec.recommender.item.RecommendedItemList;

public class Precision implements IObjectiveFunction {
	private static final String AXIS_TITLE = "Precision";
    
    @Override
    public double objectiveFunction(final ParetoObject paretoObject) {
        return objectiveFunction(paretoObject.getChromosome());
    }
    
    @Override
    public double objectiveFunction(final Chromosome chromosome)  {
        //return objectiveFunction(chromosome.getFitness());	//->Service.calculateFitness(geneticCode)
    	return objectiveFunction(chromosome.getGeneticCode());
    }														//即函数自变量入口

    @Override
    public double objectiveFunction(double geneVaue) {
        return 1;
    }

    @Override
    public double objectiveFunction(Allele[] geneticCode) {
    	double[] weights = Synthesis.weights(geneticCode);
    	double temp=0;
		try {
			AbstractRecommender recommender=RecAlgorithm.CalculatorItemList(weights);
			temp=RecAlgorithm.getPrecision(recommender)*18;
		} catch (ClassNotFoundException | IllegalAccessException | InvocationTargetException | IOException
				| LibrecException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return temp;
//    	AbstractRecommender recommender = Calculator.CalculatorItemList(weights);
//    	return Calculator.getPercision(recommender);
    }
    @Override
    public String getAxisTitle() {
        return AXIS_TITLE;
    }
}
