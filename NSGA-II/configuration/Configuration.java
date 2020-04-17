/*
 * This repository / codebase is Open Source and free for use and rewrite.
 */
package io.onclave.nsga.ii.configuration;

import io.onclave.nsga.ii.Interface.IObjectiveFunction;
import io.onclave.nsga.ii.objectivefunction.A;
import io.onclave.nsga.ii.objectivefunction.B;
import io.onclave.nsga.ii.objectivefunction.FON_1;
import io.onclave.nsga.ii.objectivefunction.FON_2;
import io.onclave.nsga.ii.objectivefunction.SCH_1;
import io.onclave.nsga.ii.objectivefunction.SCH_2;
import io.onclave.nsga.ii.objectivefunction.Precision;
import io.onclave.nsga.ii.objectivefunction.Reciprocal;
import io.onclave.nsga.ii.objectivefunction.Diversity;
import io.onclave.nsga.ii.objectivefunction.ZDT4_1;
import io.onclave.nsga.ii.objectivefunction.ZDT4_2;
import io.onclave.nsga.ii.objectivefunction.ZDT2_1;
import io.onclave.nsga.ii.objectivefunction.ZDT2_2;
import io.onclave.nsga.ii.objectivefunction.ZDT6_1;
import io.onclave.nsga.ii.objectivefunction.ZDT6_2;
import java.util.ArrayList;
import java.util.List;

/**
 * this is the Configuration file for the algorithm, where all the values are set and the initial
 * configurations are set and run.
 * to change any aspect of the algorithm, this file may be tweaked.
 * 
 * @author  Debabrata Acharya <debabrata.acharya@icloud.com>
 * @version 1.0
 * @since   0.1
 */
public class Configuration {
    
    private static final int POPULATION_SIZE = 80;
    private static final int GENERATIONS = 20;
    private static final int CHROMOSOME_LENGTH = 20;//->Synthesizes
    private static final float CROSSOVER_PROBABILITY = 0.8f;//交叉概率->Synthesizes
    private static final float MUTATION_PROBABILITY = 0.3f;//变异概率->Synthesizes0.03
    private static List<IObjectiveFunction> objectives = null;
    
    public static final double ACTUAL_MIN = 0;//->Service
    public static final double ACTUAL_MAX = Math.pow(2, CHROMOSOME_LENGTH) - 1;//->Service(setCrowdingDistance)
    public static final double NORMALIZED_MIN =0;//->Service
    public static final double NORMALIZED_MAX =1;//->Service
    public static final String DEFAULT_X_AXIS_TITLE = "x-axis";
    public static final String DEFAULT_Y_AXIS_TITLE = "y-axis";

    public static int getPOPULATION_SIZE() {
        return POPULATION_SIZE;
    }

    public static int getGENERATIONS() {
        return GENERATIONS;
    }

    public static int getCHROMOSOME_LENGTH() {
        return CHROMOSOME_LENGTH;
    }
    
    /**
     * this method sets the objective functions over which the algorithm is to operate.
     * it is a list of IObjectionFunction objects.
     */
    public static void buildObjectives() {//问题入口
        
        //List<IObjectiveFunction> newObjectives = new ArrayList<>();
    	List<IObjectiveFunction> newObjectives = new ArrayList<>();
        
        //newObjectives.add(new FON_1());//[-1/Math.sqrt(3),1/Math.sqrt(3)]
        //newObjectives.add(new FON_2());
        //newObjectives.add(new SCH_1());//[0,2]
        //newObjectives.add(new SCH_2());
       // newObjectives.add(new ZDT4_1());//[0,1]
        //newObjectives.add(new ZDT4_2());
        //newObjectives.add(new ZDT6_1());//[0,1]
        //newObjectives.add(new ZDT6_2());
        //newObjectives.add(new ZDT2_1());//[0,1]
        //newObjectives.add(new ZDT2_2());
    	newObjectives.add(new Precision());
    	//newObjectives.add(new Diversity());
    	newObjectives.add(new Reciprocal());
    	//newObjectives.add(new Test2());
    	//newObjectives.add(new A());
    	//newObjectives.add(new B());
        
        setObjectives(newObjectives);
    }

    public static List<IObjectiveFunction> getObjectives() {
        return objectives;
    }

    public static void setObjectives(List<IObjectiveFunction> objectives) {
        Configuration.objectives = objectives;
    }

    public static float getMUTATION_PROBABILITY() {
        return MUTATION_PROBABILITY;//变异概率
    }

    public static float getCROSSOVER_PROBABILITY() {
        return CROSSOVER_PROBABILITY;//交叉概率
    }
    
    public static String getXaxisTitle() {
        return getObjectives().size() > 2 ? DEFAULT_X_AXIS_TITLE : getObjectives().get(0).getAxisTitle();
    }																					//实现是在目标函数中
    
    public static String getYaxisTitle() {
        return getObjectives().size() > 2 ? DEFAULT_Y_AXIS_TITLE : getObjectives().get(1).getAxisTitle();
    }
}
