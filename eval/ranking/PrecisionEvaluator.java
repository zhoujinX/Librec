/**
 * Copyright (C) 2016 LibRec
 * <p>
 * This file is part of LibRec.
 * LibRec is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * LibRec is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with LibRec. If not, see <http://www.gnu.org/licenses/>.
 */
package net.librec.eval.ranking;

import net.librec.eval.AbstractRecommenderEvaluator;
import net.librec.math.structure.SparseMatrix;
import net.librec.recommender.item.ItemEntry;
import net.librec.recommender.item.RecommendedList;

import java.util.List;
import java.util.Set;

/**
 * PrecisionEvaluator, calculate precision@
 *
 * @author WangYuFeng
 */
public class PrecisionEvaluator extends AbstractRecommenderEvaluator {

    /**
     * Evaluate on the test set with the the list of recommended items.
     *
     * @param testMatrix
     *            the given test set
     * @param recommendedList
     *            the list of recommended items
     * @return evaluate result
     */
    public double evaluate(SparseMatrix testMatrix, RecommendedList recommendedList) {
        double totalPrecision = 0.0;
        int numUsers = testMatrix.numRows();	//用户人数
        int nonZeroNumUsers = 0;	//不为零的用户人数
        for (int userID = 0; userID < numUsers; userID++) {
            Set<Integer> testSetByUser = testMatrix.getColumnsSet(userID);//保存每个用户评过分的物品的列下标
            if (testSetByUser.size() > 0) {
            	//取出推荐列表对应用户的推荐物品列表
                List<ItemEntry<Integer, Double>> recommendListByUser = recommendedList.getItemIdxListByUserIdx(userID);
                int numHits = 0;
                int topK = this.topN <= recommendListByUser.size() ? this.topN : recommendListByUser.size();
                //取出recommendedListByUser中各个项目
                for (int indexOfItem = 0; indexOfItem < topK; indexOfItem ++) {
                    int itemID = recommendListByUser.get(indexOfItem).getKey();
                    if (testSetByUser.contains(itemID)) {	//看testSetByUser是否包含itemId
                        numHits++;
                    }
                }
                totalPrecision += numHits / (this.topN + 0.0);
                nonZeroNumUsers++;
            }
        }
        return nonZeroNumUsers > 0 ? totalPrecision / nonZeroNumUsers : 0.0d;
    }
}
