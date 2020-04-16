package net.librec.recommender;

import java.util.ArrayList;
import java.util.List;

import net.librec.common.LibrecException;
import net.librec.conf.Configuration;
import net.librec.data.model.TextDataModel;
import net.librec.eval.RecommenderEvaluator;
import net.librec.eval.ranking.DiversityEvaluator;
import net.librec.eval.ranking.PrecisionEvaluator;
import net.librec.math.algorithm.Randoms;
import net.librec.recommender.Recommender;
import net.librec.recommender.RecommenderContext;
import net.librec.recommender.cf.ranking.BPRRecommender;
import net.librec.recommender.cf.ranking.EALSRecommender;
import net.librec.recommender.cf.ranking.ListRankMFRecommender;
import net.librec.recommender.item.GenericRecommendedItem;
import net.librec.recommender.item.ItemEntry;
import net.librec.recommender.item.RecommendedItem;
import net.librec.recommender.item.RecommendedItemList;
import net.librec.similarity.PCCSimilarity;
import net.librec.similarity.RecommenderSimilarity;

//要求在主函数事先构建三个Recommender的对象
public class RecommendedCalculator {
	
	static Configuration conf = new Configuration();
	
	static RecommenderContext context = null;
	
	static RecommenderSimilarity similarity = null;
	
	static TextDataModel dataModel = null;
	
	public static Recommender recommender1 = new  BPRRecommender();
	public static Recommender recommender2 = new  EALSRecommender();
	public static Recommender recommender3 = new  ListRankMFRecommender();
	//public static RecommendedItemList itemList = 
	
	public static void recommendedCalculator() throws LibrecException {
		
		
		//建立数据模型
		conf.set("dfs.data.dir","D:/Librec/librec/data");//设置读取数据路径
		Randoms.seed(1);
		dataModel = new TextDataModel(conf);
		dataModel.buildDataModel();//建立数据模型
						
		//建立推荐系统环境
		context = new RecommenderContext(conf,dataModel);
						
		//创建相似度矩阵
		conf.set("rec.recommender.similarity.key", "item");//以item为维度建立相似矩阵
		similarity = new PCCSimilarity();
		similarity = new PCCSimilarity();
		similarity.buildSimilarityMatrix(dataModel);
		context.setSimilarity(similarity);
		
		//推荐算法运行
		recommender1.setContext(context);
		recommender1.recommend(context);
		recommender2.setContext(context);
		recommender2.recommend(context);
		recommender3.setContext(context);
		recommender3.recommend(context);
	}
	
	//通过权重计算出对应的推荐列表
	public static RecommendedItemList CalculatorItemList(double[] weighted,Recommender recommender1,Recommender recommender2,Recommender recommender3) {
		//将要用到的参数
		List<RecommendedItem> recommendedItemList1 = null;
		List<RecommendedItem> recommendedItemList2 = null;
		List<RecommendedItem> recommendedItemList3 = null;
		
		//综合推荐结果
		recommendedItemList1 = recommender1.getRecommendedList();
		recommendedItemList2 = recommender2.getRecommendedList();
		recommendedItemList3 = recommender3.getRecommendedList();
		
		List<RecommendedItem> newRcommendedItemList = new ArrayList<>();
		GenericRecommendedItem tempRecommendedItem = null;
		
		for(int i = 0;i < recommendedItemList1.size();i ++) {
			tempRecommendedItem = new GenericRecommendedItem(
					recommendedItemList1.get(i).getUserId(),
					recommendedItemList1.get(i).getItemId(),
						(recommendedItemList1.get(i).getValue() * weighted[0]) +
						(recommendedItemList2.get(i).getValue() * weighted[1]) +
						(recommendedItemList3.get(i).getValue() * weighted[2]));
			newRcommendedItemList.add(tempRecommendedItem);
		}
		
		//构建RecommendedItemList对象并返回
		RecommendedItemList recommendedItemList = new RecommendedItemList(recommendedItemList1.size());
		
		ArrayList<ItemEntry<Integer,Double>> itemList = new ArrayList<ItemEntry<Integer,Double>>();
		for(int i = 1;i <= recommendedItemList1.size();i ++) {
			for(RecommendedItem recommendedItem : newRcommendedItemList)
			{
				if(Integer.parseInt(recommendedItem.getUserId()) == i) {
					ItemEntry<Integer,Double> item = new ItemEntry<Integer,Double>(Integer.parseInt(recommendedItem.getItemId()),
							recommendedItem.getValue());
					itemList.add(item);
				}
			}
			recommendedItemList.addItemIdxList(i, itemList);
			itemList.clear();
		}
		
		return recommendedItemList;
	}
	
	public static double[] Evaluator(RecommendedItemList recommendedItemList) {
		//用于存储P、D值
		double []eval_r=new double[2];
		
		//System.out.println(recommendedItemList.getEntryValue(196, 382));
		
		//计算P、D值
		conf.set("rec.recommender.isranking","true");
	    conf.set("rec.eval.enable","true");
	    RecommenderEvaluator percisionEvaluator = new PrecisionEvaluator();
		percisionEvaluator.setTopN(10);
		double precision = percisionEvaluator.evaluate(context, recommendedItemList);
		//System.out.println(precision);
		eval_r[0]=precision;
		
		/*String[] similarityKeys = conf.getStrings("rec.recommender.similarities");//item键
	        if (similarityKeys != null && similarityKeys.length > 0) {
	            for(int i = 0; i< similarityKeys.length; i++) {
	                    conf.set("rec.recommender.similarity.key", similarityKeys[i]);
	                    similarity.buildSimilarityMatrix(dataModel);
	                    if(i == 0){
	                        context.setSimilarity(similarity);
	                    }
	                    context.addSimilarities(similarityKeys[i], similarity);
	        }
	    }
	
		RecommenderEvaluator diversityEvaluator =new DiversityEvaluator();
		diversityEvaluator.setTopN(10);
		double diversity = diversityEvaluator.evaluate(context, recommendedItemList);
		System.out.println(diversity);
		eval_r[1]=diversity;*/
		
		return eval_r;
	}
}
