package Demo;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.beanutils.BeanUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import net.librec.common.LibrecException;
import net.librec.conf.Configuration;
import net.librec.eval.RecommenderEvaluator;
import net.librec.eval.ranking.DiversityEvaluator;
import net.librec.eval.ranking.PrecisionEvaluator;
import net.librec.eval.ranking.RecallEvaluator;
import net.librec.job.RecommenderJob;
import net.librec.recommender.AbstractRecommender;
import net.librec.recommender.Recommender;
import net.librec.recommender.RecommenderContext;
import net.librec.recommender.item.GenericRecommendedItem;
import net.librec.recommender.item.ItemEntry;
import net.librec.recommender.item.RecommendedItem;
import net.librec.recommender.item.RecommendedItemList;
import net.librec.recommender.item.RecommendedList;
import net.librec.util.ReflectionUtil;

public class RecAlgorithm {

	private static AbstractRecommender recommender1;
	private static AbstractRecommender recommender2;
	private static AbstractRecommender recommender;
	
	private static Configuration conf1;
	private static Configuration conf2;
	private Configuration conf;
	
	 static RecommenderJob job1;
	 static RecommenderJob job2;
	static RecommenderJob job;
	
	static RecommenderContext context1 = null;
	static RecommenderContext context2 = null;
	
	public static void RecBuildModel()  throws ClassNotFoundException, LibrecException, IOException, IllegalAccessException, InvocationTargetException {
		
		conf1=new Configuration();
		conf2=new Configuration();
		
        Configuration.Resource resource1 = new Configuration.Resource("rec/cf/ranking/ranksgd-test.properties");
        conf1.addResource(resource1);
        conf1.setInt("rec.recommender.ranking.topn", 10);
        job1 = new RecommenderJob(conf1);
        job1.generateDataModel();
        context1 = new RecommenderContext(conf1, job1.getDataModel());
        job1.generateSimilarity(context1);
        recommender1 = (AbstractRecommender) ReflectionUtil.newInstance((Class<Recommender>) job1.getRecommenderClass(), conf1);
        recommender1.recommend(context1);
        
//        job1.executeEvaluator(recommender1);
//        List<RecommendedItem> recommendedList1 = recommender1.getRecommendedList();
//        //recommendedList = filterResult(recommendedList);
////        job1.saveResult(recommendedList1);

        Configuration.Resource resource2 = new Configuration.Resource("rec/cf/ranking/gbpr-test.properties");
        conf2.addResource(resource2);
        conf2.setInt("rec.recommender.ranking.topn", 10);
        job2 = new RecommenderJob(conf2);
        job2.generateDataModel();
        context2 = new RecommenderContext(conf2, job2.getDataModel());
        job2.generateSimilarity(context2);
        recommender2 = (AbstractRecommender) ReflectionUtil.newInstance((Class<Recommender>) job2.getRecommenderClass(), conf2);
        recommender2.recommend(context2);
        
//        job2.executeEvaluator(recommender2);
//        List<RecommendedItem> recommendedList2 = recommender2.getRecommendedList();
////        recommendedList = filterResult(recommendedList);
////        job2.saveResult(recommendedList2);
//        
//        List<RecommendedItem> newRecommendedItemList = new ArrayList<>();//加权后的的推荐列表合格
//		  GenericRecommendedItem tempRecommendedItem = null;
//	     
//        for(int i=0;i<recommendedList1.size();i++) {
//	    	// RecommendedItem tempItem;
//	    	 double newValue=recommendedList1.get(i).getValue()*0.2+recommendedList2.get(i).getValue()*0.8;
//	    	 int temUserId=Integer.parseInt(recommendedList1.get(i).getUserId());
//	    	 int temItemId=Integer.parseInt(recommendedList1.get(i).getItemId());
//	    	 //recommendedItemList.addUserItemIdx(temUserId, temItemId, newValue);
//	    	 tempRecommendedItem = new GenericRecommendedItem(recommendedList1.get(i).getUserId(),recommendedList1.get(i).getItemId(),newValue);//构造函数
//			 newRecommendedItemList.add(tempRecommendedItem);
//	     }
//
//        RecommendedItemList recommendedItemList=new RecommendedItemList(newRecommendedItemList.size(),newRecommendedItemList.size()+1);
//		  for (RecommendedItem recommendedItem : newRecommendedItemList) { 
//			  int userId=Integer.parseInt(recommendedItem.getUserId()); 
//			  int itemId=Integer.parseInt(recommendedItem.getItemId());
//
//			  recommendedItemList.addUserItemIdx(userId-1, itemId,recommendedItem.getValue()); 
//		  }
//		  
//		  for (int userIdx = 0; userIdx < 943; ++userIdx) {
//		    	 recommender1.recommendedList.topNRankItemsByUser(userIdx, 10);
//			    }
//		    job1.saveResult(recommender1.getRecommendedList());
//
//		    for (int userIdx = 0; userIdx < 943; ++userIdx) {
//		    	 recommender2.recommendedList.topNRankItemsByUser(userIdx, 10);
//			    }
//		    job2.saveResult(recommender2.getRecommendedList());
////      recommender=recommender1;
////		recommender = (AbstractRecommender) ReflectionUtil.newInstance((Class<Recommender>) job1.getRecommenderClass(), conf1);
////		recommender.recommendedList=new RecommendedItemList(942,943); 
////		BeanUtils.copyProperties(recommender1.recommendedList,recommendedItemList);
//		
//		recommender1.recommendedList=recommendedItemList;
////		recommender1.setRecommendedList(recommendedItemList);
////		recommender1.recommendedList=(RecommendedItemList) BeanUtils.cloneBean(recommendedItemList);
//		for (int userIdx = 0; userIdx < 943; ++userIdx) {
//	    	 recommender1.recommendedList.topNRankItemsByUser(userIdx, 10);
//		    }
//	    job2.saveResult2(recommender1.getRecommendedList());
    }
	
	public static AbstractRecommender CalculatorItemList(double[] weighted) throws ClassNotFoundException, IOException, LibrecException, IllegalAccessException, InvocationTargetException {
		 
		 if(weighted[0]==1.0) return recommender1;
		 if(weighted[1]==1.0) return recommender2;
		
		 List<RecommendedItem> recommendedList1 = recommender1.getRecommendedList();
	     recommendedList1 = job1.filterResult(recommendedList1);
	     List<RecommendedItem> recommendedList2 = recommender2.getRecommendedList();
	     recommendedList2 = job2.filterResult(recommendedList2);
	     //List<RecommendedItem> recommendedList=null;
	     
	     List<RecommendedItem> newRecommendedItemList = new ArrayList<>();//加权后的的推荐列表合格
		 GenericRecommendedItem tempRecommendedItem = null;
		     
	        for(int i=0;i<recommendedList1.size();i++) {
		    	// RecommendedItem tempItem;
		    	 double newValue=recommendedList1.get(i).getValue()*weighted[0]+recommendedList2.get(i).getValue()*weighted[1];
		    	 int temUserId=Integer.parseInt(recommendedList1.get(i).getUserId());
		    	 int temItemId=Integer.parseInt(recommendedList1.get(i).getItemId());
		    	 //recommendedItemList.addUserItemIdx(temUserId, temItemId, newValue);
		    	 tempRecommendedItem = new GenericRecommendedItem(recommendedList1.get(i).getUserId(),recommendedList1.get(i).getItemId(),newValue);//构造函数
				 newRecommendedItemList.add(tempRecommendedItem);
		     }
	        
//	        -------------下面addUserItemIdx有问题，上面的加权综合推荐列表无误
	          RecommendedItemList recommendedItemList=new RecommendedItemList(942,943);
			  for (RecommendedItem recommendedItem : newRecommendedItemList) { 
				  int userId=Integer.parseInt(recommendedItem.getUserId()); 
				  int itemId=Integer.parseInt(recommendedItem.getItemId());
				  
			  	  recommendedItemList.addUserItemIdx(userId-1, itemId,recommendedItem.getValue()); 
			  }
			 
	        recommender=recommender1;
//	        BeanUtils.copyProperties(recommender.recommendedList, recommendedItemList);
		    recommender.recommendedList=recommendedItemList;
//	        recommender.recommendedList=(RecommendedList) BeanUtils.cloneBean(recommendedItemList);
		    for (int userIdx = 0; userIdx < 943; ++userIdx) {
	    	 recommender.recommendedList.topNRankItemsByUser(userIdx, 10);
		    }
//	     job = new RecommenderJob(recommender.conf);
//	     job1.saveResult2(recommender.getRecommendedList());
	     //writeFile(recommender);
	    return recommender;
	    

	}
	
	protected final Log LOG = LogFactory.getLog(RecAlgorithm.class);
	
public static double getPrecision2(RecommendedItemList recommendedItemList) throws ClassNotFoundException, IOException, LibrecException{
		
	conf1.set("rec.recommender.isranking","true");
    conf1.set("rec.eval.enable","true");
    RecommenderEvaluator percisionEvaluator = new PrecisionEvaluator();
	percisionEvaluator.setTopN(10);
	double precision = percisionEvaluator.evaluate(context2, recommendedItemList);
	return precision;
	}

public static double getRecall2(RecommendedItemList recommendedItemList) throws ClassNotFoundException, IOException, LibrecException{
	
	RecommenderEvaluator diversityEvaluator =new RecallEvaluator();
	diversityEvaluator.setTopN(10);
	double recall = diversityEvaluator.evaluate(context2, recommendedItemList);
	return recall;
	}
	
	public static double getPrecision(AbstractRecommender recommender) throws ClassNotFoundException, IOException, LibrecException{
		
		job = new RecommenderJob(conf1);
		RecommenderEvaluator evaluator = (RecommenderEvaluator) ReflectionUtil.newInstance(job.getEvaluatorClass("precision"), null);
        evaluator.setTopN(recommender.conf.getInt("rec.recommender.ranking.topn", 10));
        double evalValue = recommender.evaluate(evaluator);
//        job.LOG.info("Evaluator info:" + evaluator.getClass().getSimpleName() + " is " + evalValue);
		return evalValue;
	}
	
	public static double getRecall(AbstractRecommender recommender) throws ClassNotFoundException, IOException, LibrecException{
		job = new RecommenderJob(conf1);
	RecommenderEvaluator evaluator = (RecommenderEvaluator) ReflectionUtil.newInstance(job.getEvaluatorClass("novelty"), null);
        evaluator.setTopN(recommender.conf.getInt("rec.recommender.ranking.topn", 10));
        double evalValue = recommender.evaluate(evaluator);
//        job.LOG.info("Evaluator info:" + evaluator.getClass().getSimpleName() + " is " + evalValue);
		return evalValue;
	}
	
	public static void writeFile(Recommender recommender) throws IOException {
		List<RecommendedItem> recommendedList = recommender.getRecommendedList();
		FileWriter fw = new FileWriter("D:\\Librec\\librec\\result\\movielens\\ml-100k\\result.txt",true);	
    	
    	for(RecommendedItem item:recommendedList) {
    		
    		fw.write(item.getUserId()+"\t" +item.getItemId()+"\t"+item.getValue()+"\r\n");
    	}
    	fw.write("\r\n");
    	fw.close();
    }
	
	public static void main(String[] args) throws ClassNotFoundException, LibrecException, IOException, IllegalAccessException, InvocationTargetException, NoSuchMethodException, InstantiationException {
		RecAlgorithm r=new RecAlgorithm();
		r.RecBuildModel();
		double[] weighted= {0.5,0.5};
		AbstractRecommender recommender=r.CalculatorItemList(weighted);
		System.out.println(r.getPrecision(recommender)+"\t"+r.getRecall(recommender));
	}

}
