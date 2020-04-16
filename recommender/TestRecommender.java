package net.librec.recommender;
import java.util.ArrayList;
import java.util.List;
import net.librec.common.LibrecException;
import net.librec.conf.Configuration;
import net.librec.data.model.TextDataModel;
import net.librec.data.splitter.KCVDataSplitter;
import net.librec.data.splitter.LOOCVDataSplitter;
import net.librec.eval.RecommenderEvaluator;
import net.librec.eval.ranking.DiversityEvaluator;
import net.librec.eval.ranking.PrecisionEvaluator;
import net.librec.math.algorithm.Randoms;
import net.librec.recommender.Recommender;
import net.librec.recommender.RecommenderContext;
import net.librec.recommender.cf.ranking.BPRRecommender;
import net.librec.recommender.cf.ranking.ListRankMFRecommender;
import net.librec.recommender.cf.ranking.RankSGDRecommender;
import net.librec.recommender.item.GenericRecommendedItem;
import net.librec.recommender.item.ItemEntry;
import net.librec.recommender.item.RecommendedItem;
import net.librec.recommender.item.RecommendedItemList;
import net.librec.recommender.item.RecommendedList;
import net.librec.similarity.PCCSimilarity;
import net.librec.similarity.RecommenderSimilarity;
import net.librec.util.FileUtil;

public class TestRecommender {
	
	static Configuration conf = new Configuration();
	static RecommenderContext context = null;
	public static TextDataModel dataModel1 = null;
	
	public static Recommender[] recommenders = new Recommender[3];////////
	
	/*public static void main(String[] args) throws Exception
	{	
		double[] weighteds = {0.2,0.3,0.5};//0.2,0.3,0.5
		TextDataModel dataModel = goujian();
		recommender(weighteds,dataModel);
		//Initiate();
	}*/
	
	static public void/*TextDataModel*/ goujian() throws LibrecException{
		//建立数据模型
		conf.set("dfs.data.dir","D:/Librec/librec/data");//设置读取数据路径../data
		conf.set("data.input.path", "movielens/ml-100k/ratings.txt");//数据集设置为ml-100k
		Randoms.seed(1);
		TextDataModel dataModel = new TextDataModel(conf);
		dataModel.buildDataModel();//建立数据模型

		//建立推荐系统环境
		context = new RecommenderContext(conf,dataModel);
       
		//创建相似度矩阵
		conf.set("rec.recommender.similarity.key", "item");//以item为维度建立相似矩阵
		conf.set("rec.recommender.similarities","item");
		RecommenderSimilarity similarity = new PCCSimilarity();


		 String[] similarityKeys = conf.getStrings("rec.recommender.similarities");//item键
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

		//评估器配置参数
		conf.set("rec.recommender.isranking","true");
        conf.set("rec.eval.enable","true");
        
        //return dataModel;
	}
	
	public static void ThreeRecommender(TextDataModel dataModel) throws LibrecException {
		//推荐算法运行
				//Recommender[] recommenders = new Recommender[3];
				recommenders[0] = new BPRRecommender();
				recommenders[1] = new ListRankMFRecommender();
				recommenders[2] = new RankSGDRecommender();
				for(int i = 0;i <3;i ++) {
					recommenders[i].setContext(context);
					recommenders[i].recommend(context);	
				}

	}
	
	public static double[] recommender(double[] weighteds,Recommender[] recommenders,TextDataModel dataModel) throws LibrecException//double[] weighteds,TextDataModel dataModel
	{	
		double[]  PD= new double[2];
		/*//推荐算法运行
		Recommender[] recommenders = new Recommender[3];
		recommenders[0] = new BPRRecommender();
		recommenders[1] = new ListRankMFRecommender();
		recommenders[2] = new RankSGDRecommender();
		for(int i = 0;i <3;i ++) {
			recommenders[i].setContext(context);
			recommenders[i].recommend(context);	
		}
		*/		
		//初始化权重
		double[] weighted = weighteds;
				
		//设置输出路径
		List<RecommendedItem> BPRRecommendedItemList = null;
		List<RecommendedItem> ListRankMFRecommendedItemList = null;
		List<RecommendedItem> RankSGDRecommendedItemList = null;
		StringBuilder sb = new StringBuilder();
		String resultData = null;
		String algoSimpleName = "result";
		String outputPath = conf.get("dfs.result.dir") + "/" + conf.get("data.input.path") + "-" + algoSimpleName + "-output/" + algoSimpleName;//输出路径
			      if (null != dataModel && (dataModel.getDataSplitter() instanceof KCVDataSplitter || dataModel.getDataSplitter() instanceof LOOCVDataSplitter) && null != conf.getInt("data.splitter.cv.index")) {
			          outputPath = outputPath + "-" + String.valueOf(conf.getInt("data.splitter.cv.index"));
			      }
			      //System.out.println("Result path is " + outputPath);

		//综合三个推荐结果     
		BPRRecommendedItemList = recommenders[0].getRecommendedList();
		ListRankMFRecommendedItemList = recommenders[1].getRecommendedList();
		RankSGDRecommendedItemList = recommenders[2].getRecommendedList();
		
		List<RecommendedItem> newRcommendedItemList = new ArrayList<>();
		GenericRecommendedItem tempRecommendedItem = null;
		
		for(int i = 0;i < BPRRecommendedItemList.size();i++) {
			tempRecommendedItem = new GenericRecommendedItem(
					BPRRecommendedItemList.get(i).getUserId(),
					BPRRecommendedItemList.get(i).getItemId(),
						(BPRRecommendedItemList.get(i).getValue() * weighted[0]) +
						(ListRankMFRecommendedItemList.get(i).getValue() * weighted[1])+
						(RankSGDRecommendedItemList.get(i).getValue() * weighted[2]));
			newRcommendedItemList.add(tempRecommendedItem);
		}
		
		//评估
		RecommendedItemList recommendedItemList = new RecommendedItemList(BPRRecommendedItemList.size());
		ArrayList<ItemEntry<Integer,Double>> itemList = new ArrayList<ItemEntry<Integer,Double>>();
		for(int i = 1;i <= BPRRecommendedItemList.size();i++) {
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
		
		/*for(RecommendedItem recommendedItem : newRcommendedItemList)
		{
			recommendedItemList.addUserItemIdx(Integer.parseInt(recommendedItem.getUserId()), 
					Integer.parseInt(recommendedItem.getItemId()), 
					recommendedItem.getValue());	
		}*/
		
		//结果分析P、D值
		RecommenderEvaluator percisionEvaluator = new PrecisionEvaluator();
		percisionEvaluator.setTopN(10);//设置TopN,可变
		double precision = percisionEvaluator.evaluate(context, recommendedItemList);
		sb.append(precision).append("\t");//sb.append("Precision:").append(precision).append("\n");
		System.out.println(precision);
		PD[0] = precision;
	
		RecommenderEvaluator diversityEvaluator =new DiversityEvaluator();
		diversityEvaluator.setTopN(10);
		double diversity = diversityEvaluator.evaluate(context, recommendedItemList);
		sb.append(diversity).append("\n");//sb.append("Diversity:").append(diversity).append("\n");
		System.out.println(diversity);
		PD[0] = diversity;
		
		//保存结果进入文件
		resultData = sb.toString();
		try {
				FileUtil.writeString(outputPath, resultData);//保存输出文件
			} catch (Exception e) {
				e.printStackTrace();
		}
		return PD;
	}
	

}
