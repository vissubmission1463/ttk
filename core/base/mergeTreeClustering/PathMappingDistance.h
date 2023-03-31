/// \ingroup base
/// \class PathMappingDistance
/// \author Florian Wetzels (wetzels@cs.uni-kl.de)
/// \date 2022.
///
/// This module defines the %PathMappingDistance class that computes distances
/// between two merge trees.
///
/// \b Related \b publication \n
/// "A Deformation-based Edit Distance for Merge Trees" \n
/// Florian Wetzels, Christoph Garth. \n
/// TopoInVis 2022.

#pragma once

#include <set>
#include <vector>

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <tuple>
#include <vector>

// ttk common includes
#include "MergeTreeBase.h"
#include <AssignmentAuction.h>
#include <AssignmentExhaustive.h>
#include <AssignmentMunkres.h>
#include <Debug.h>
#include <FTMTree_MT.h>

namespace ttk {

  class PathMappingDistance : virtual public Debug, public MergeTreeBase {

  private:
    int baseMetric_ = 0;
    int assignmentSolverID_ = 0;
    bool squared_ = false;
    bool computeMapping_ = false;

    bool preprocess_ = true;

    template <class dataType>
    inline dataType editCost_Persistence(int n1,
                                         int p1,
                                         int n2,
                                         int p2,
                                         ftm::FTMTree_MT *tree1,
                                         ftm::FTMTree_MT *tree2) {
      dataType d;
      if(n1 < 0) {
        dataType b1 = tree2->getValue<dataType>(n2);
        dataType d1 = tree2->getValue<dataType>(p2);
        d = (d1 > b1) ? (d1 - b1) : (b1 - d1);
      } else if(n2 < 0) {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        d = (d1 > b1) ? (d1 - b1) : (b1 - d1);
      } else {
        dataType b1 = tree1->getValue<dataType>(n1);
        dataType d1 = tree1->getValue<dataType>(p1);
        dataType b2 = tree2->getValue<dataType>(n2);
        dataType d2 = tree2->getValue<dataType>(p2);
        dataType dist1 = (d1 > b1) ? (d1 - b1) : (b1 - d1);
        dataType dist2 = (d2 > b2) ? (d2 - b2) : (b2 - d2);
        d = (dist1 > dist2) ? (dist1 - dist2) : (dist2 - dist1);
      }
      return squared_ ? d * d : d;
    }

    template <class dataType>
    void traceMapping_path(
      ftm::FTMTree_MT *tree1,
      ftm::FTMTree_MT *tree2,
      int curr1,
      int l1,
      int curr2,
      int l2,
      std::vector<std::vector<int>> &predecessors1,
      std::vector<std::vector<int>> &predecessors2,
      int depth1,
      int depth2,
      std::vector<dataType> &memT,
      std::vector<std::pair<std::pair<ftm::idNode, ftm::idNode>,
                            std::pair<ftm::idNode, ftm::idNode>>> &mapping) {

      //===============================================================================
      // If both trees not empty, find optimal edit operation
      std::vector<ftm::idNode> children1;
      tree1->getChildren(curr1, children1);
      std::vector<ftm::idNode> children2;
      tree2->getChildren(curr2, children2);
      int parent1 = predecessors1[curr1][predecessors1[curr1].size() - l1];
      int parent2 = predecessors2[curr2][predecessors2[curr2].size() - l2];

      size_t nn1 = tree1->getNumberOfNodes();
      size_t nn2 = tree2->getNumberOfNodes();
      size_t dim1 = 1;
      size_t dim2 = (nn1 + 1) * dim1;
      size_t dim3 = (depth1 + 1) * dim2;
      size_t dim4 = (nn2 + 1) * dim3;

      //---------------------------------------------------------------------------
      // If both trees only have one branch, return edit cost between
      // the two branches
      if(tree1->getNumberOfChildren(curr1) == 0
         and tree2->getNumberOfChildren(curr2) == 0) {
        mapping.emplace_back(std::make_pair(
          std::make_pair(curr1, parent1), std::make_pair(curr2, parent2)));
        return;
      }
      //---------------------------------------------------------------------------
      // If first tree only has one branch, try all decompositions of
      // second tree
      else if(tree1->getNumberOfChildren(curr1) == 0) {
        for(auto child2_mb : children2) {
          dataType d_
            = memT[curr1 + l1 * dim2 + child2_mb * dim3 + (l2 + 1) * dim4];
          for(auto child2 : children2) {
            if(child2 == child2_mb) {
              continue;
            }
            d_ += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
          }
          if(memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] == d_) {
            traceMapping_path(tree1, tree2, curr1, l1, child2_mb, l2 + 1,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            return;
          }
        }
      }
      //---------------------------------------------------------------------------
      // If second tree only has one branch, try all decompositions of
      // first tree
      else if(tree2->getNumberOfChildren(curr2) == 0) {
        dataType d = std::numeric_limits<dataType>::max();
        for(auto child1_mb : children1) {
          dataType d_
            = memT[child1_mb + (l1 + 1) * dim2 + curr2 * dim3 + l2 * dim4];
          for(auto child1 : children1) {
            if(child1 == child1_mb) {
              continue;
            }
            d_ += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
          }
          d = std::min(d, d_);
          if(memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] == d_) {
            traceMapping_path(tree1, tree2, child1_mb, l1 + 1, curr2, l2,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            return;
          }
        }
      }
      //---------------------------------------------------------------------------
      // If both trees have more than one branch, try all decompositions
      // of both trees
      else {
        //-----------------------------------------------------------------------
        // Try all possible main branches of first tree (child1_mb) and
        // all possible main branches of second tree (child2_mb) Then
        // try all possible matchings of subtrees
        if(tree1->getNumberOfChildren(curr1) == 2
           && tree2->getNumberOfChildren(curr2) == 2) {
          int child11 = children1[0];
          int child12 = children1[1];
          int child21 = children2[0];
          int child22 = children2[1];
          if(memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4]
             == memT[child11 + 1 * dim2 + child21 * dim3 + 1 * dim4]
                  + memT[child12 + 1 * dim2 + child22 * dim3 + 1 * dim4]
                  + editCost_Persistence<dataType>(
                    curr1, parent1, curr2, parent2, tree1, tree2)) {
            mapping.emplace_back(std::make_pair(
              std::make_pair(curr1, parent1), std::make_pair(curr2, parent2)));
            traceMapping_path(tree1, tree2, child11, 1, child21, 1,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            traceMapping_path(tree1, tree2, child12, 1, child22, 1,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            return;
          }
          if(memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4]
             == memT[child11 + 1 * dim2 + child22 * dim3 + 1 * dim4]
                  + memT[child12 + 1 * dim2 + child21 * dim3 + 1 * dim4]
                  + editCost_Persistence<dataType>(
                    curr1, parent1, curr2, parent2, tree1, tree2)) {
            mapping.emplace_back(std::make_pair(
              std::make_pair(curr1, parent1), std::make_pair(curr2, parent2)));
            traceMapping_path(tree1, tree2, child11, 1, child22, 1,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            traceMapping_path(tree1, tree2, child12, 1, child21, 1,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            return;
          }
        } else {
          auto f = [&](int r, int c) {
            size_t c1
              = r < tree1->getNumberOfChildren(curr1) ? children1[r] : nn1;
            size_t c2
              = c < tree2->getNumberOfChildren(curr2) ? children2[c] : nn2;
            int l1_ = c1 == nn1 ? 0 : 1;
            int l2_ = c2 == nn2 ? 0 : 1;
            return memT[c1 + l1_ * dim2 + c2 * dim3 + l2_ * dim4];
          };
          int size = std::max(tree1->getNumberOfChildren(curr1),
                              tree2->getNumberOfChildren(curr2))
                     + 1;
          auto costMatrix = std::vector<std::vector<dataType>>(
            size, std::vector<dataType>(size, 0));
          std::vector<MatchingType> matching;
          for(int r = 0; r < size; r++) {
            for(int c = 0; c < size; c++) {
              costMatrix[r][c] = f(r, c);
            }
          }

          AssignmentSolver<dataType> *assignmentSolver;
          AssignmentExhaustive<dataType> solverExhaustive;
          AssignmentMunkres<dataType> solverMunkres;
          AssignmentAuction<dataType> solverAuction;
          switch(assignmentSolverID_) {
            case 1:
              solverExhaustive = AssignmentExhaustive<dataType>();
              assignmentSolver = &solverExhaustive;
              break;
            case 2:
              solverMunkres = AssignmentMunkres<dataType>();
              assignmentSolver = &solverMunkres;
              break;
            case 0:
            default:
              solverAuction = AssignmentAuction<dataType>();
              assignmentSolver = &solverAuction;
          }
          assignmentSolver->setInput(costMatrix);
          assignmentSolver->setBalanced(true);
          assignmentSolver->run(matching);
          dataType d_ = editCost_Persistence<dataType>(
            curr1, parent1, curr2, parent2, tree1, tree2);
          for(auto m : matching)
            d_ += std::get<2>(m);
          if(memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] == d_) {
            mapping.emplace_back(std::make_pair(
              std::make_pair(curr1, parent1), std::make_pair(curr2, parent2)));
            for(auto m : matching) {
              int n1 = std::get<0>(m) < tree1->getNumberOfChildren(curr1)
                         ? children1[std::get<0>(m)]
                         : -1;
              int n2 = std::get<1>(m) < tree2->getNumberOfChildren(curr2)
                         ? children2[std::get<1>(m)]
                         : -1;
              if(n1 >= 0 && n2 >= 0)
                traceMapping_path(tree1, tree2, n1, 1, n2, 1, predecessors1,
                                  predecessors2, depth1, depth2, memT, mapping);
            }
            return;
          }
        }
        //-----------------------------------------------------------------------
        // Try to continue main branch on one child of first tree and
        // delete all other subtrees Then match continued branch to
        // current branch in second tree
        for(auto child1_mb : children1) {
          dataType d_
            = memT[child1_mb + (l1 + 1) * dim2 + curr2 * dim3 + l2 * dim4];
          for(auto child1 : children1) {
            if(child1 == child1_mb) {
              continue;
            }
            d_ += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
          }
          if(memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] == d_) {
            traceMapping_path(tree1, tree2, child1_mb, l1 + 1, curr2, l2,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            return;
          }
        }
        //-----------------------------------------------------------------------
        // Try to continue main branch on one child of second tree and
        // delete all other subtrees Then match continued branch to
        // current branch in first tree
        for(auto child2_mb : children2) {
          dataType d_
            = memT[curr1 + l1 * dim2 + child2_mb * dim3 + (l2 + 1) * dim4];
          for(auto child2 : children2) {
            if(child2 == child2_mb) {
              continue;
            }
            d_ += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
          }
          if(memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] == d_) {
            traceMapping_path(tree1, tree2, curr1, l1, child2_mb, l2 + 1,
                              predecessors1, predecessors2, depth1, depth2,
                              memT, mapping);
            return;
          }
        }
      }
    }

  public:
    PathMappingDistance() {
      this->setDebugMsgPrefix(
        "MergeTreeDistance"); // inherited from Debug: prefix will be printed at
                              // the beginning of every msg
    }
    ~PathMappingDistance() override = default;

    void setBaseMetric(int m) {
      baseMetric_ = m;
    }

    void setAssignmentSolver(int assignmentSolver) {
      assignmentSolverID_ = assignmentSolver;
    }

    void setSquared(bool s) {
      squared_ = s;
    }

    void setComputeMapping(bool m) {
      computeMapping_ = m;
    }

    void setPreprocess(bool p) {
      preprocess_ = p;
    }

    template <class dataType>
    dataType editDistance_path(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, std::vector<std::pair<std::pair<ftm::idNode, ftm::idNode>,std::pair<ftm::idNode, ftm::idNode>>> *outputMatching) {

      // optional preprocessing

      if(preprocess_){
        preprocessTree<dataType>(tree1, true);
        preprocessTree<dataType>(tree2, true);

        // - Delete null persistence pairs and persistence thresholding
        persistenceThresholding<dataType>(tree1, persistenceThreshold_);
        persistenceThresholding<dataType>(tree2, persistenceThreshold_);

        // - Merge saddle points according epsilon
        if(not isPersistenceDiagram_) {
          treesNodeCorr_.resize(2);
          if(epsilonTree1_ != 0){
            std::vector<std::vector<ftm::idNode>> treeNodeMerged1( tree1->getNumberOfNodes() );
            mergeSaddle<dataType>(tree1, epsilonTree1_, treeNodeMerged1);
            for(unsigned int i=0; i<treeNodeMerged1.size(); i++){
              for(auto j : treeNodeMerged1[i]){
                auto nodeToDelete = tree1->getNode(j)->getOrigin();
                tree1->getNode(j)->setOrigin(i);
                tree1->getNode(nodeToDelete)->setOrigin(-1);
              }
            }
            ftm::cleanMergeTree<dataType>(tree1, treesNodeCorr_[0], true);
          }
          else{
            std::vector<ttk::SimplexId> nodeCorr1(tree1->getNumberOfNodes());
            for(unsigned int i=0; i<nodeCorr1.size(); i++) nodeCorr1[i] = i;
            treesNodeCorr_[0] = nodeCorr1;
          }
          if(epsilonTree2_ != 0){
            std::vector<std::vector<ftm::idNode>> treeNodeMerged2( tree2->getNumberOfNodes() );
            mergeSaddle<dataType>(tree2, epsilonTree2_, treeNodeMerged2);
            for(unsigned int i=0; i<treeNodeMerged2.size(); i++){
              for(auto j : treeNodeMerged2[i]){
                auto nodeToDelete = tree2->getNode(j)->getOrigin();
                tree2->getNode(j)->setOrigin(i);
                tree2->getNode(nodeToDelete)->setOrigin(-1);
              }
            }
            ftm::cleanMergeTree<dataType>(tree2, treesNodeCorr_[1], true);
          }
          else{
            std::vector<ttk::SimplexId> nodeCorr2(tree2->getNumberOfNodes());
            for(unsigned int i=0; i<nodeCorr2.size(); i++) nodeCorr2[i] = i;
            treesNodeCorr_[1] = nodeCorr2;
          }
        }
        
        if(deleteMultiPersPairs_)
          deleteMultiPersPairs<dataType>(tree1, false);
        if(deleteMultiPersPairs_)
          deleteMultiPersPairs<dataType>(tree2, false);
      }
      
      // compute preorder of both trees (necessary for bottom-up dynamic programming)

      std::vector<std::vector<int>> predecessors1(tree1->getNumberOfNodes());
      std::vector<std::vector<int>> predecessors2(tree2->getNumberOfNodes());
      int rootID1 = tree1->getRoot();
      int rootID2 = tree2->getRoot();
      std::vector<int> preorder1(tree1->getNumberOfNodes());
      std::vector<int> preorder2(tree2->getNumberOfNodes());

      int depth1 = 0;
      int depth2 = 0;
      std::stack<int> stack;
      stack.push(rootID1);
      int count = tree1->getNumberOfNodes() - 1;
      while(!stack.empty()) {
        int nIdx = stack.top();
        stack.pop();
        preorder1[count] = nIdx;
        count--;
        depth1 = std::max((int)predecessors1[nIdx].size(), depth1);
        std::vector<ftm::idNode> children;
        tree1->getChildren(nIdx, children);
        for(int cIdx : children) {
          stack.push(cIdx);
          predecessors1[cIdx].reserve(predecessors1[nIdx].size() + 1);
          predecessors1[cIdx].insert(predecessors1[cIdx].end(),
                                     predecessors1[nIdx].begin(),
                                     predecessors1[nIdx].end());
          predecessors1[cIdx].push_back(nIdx);
        }
      }
      stack.push(rootID2);
      count = tree2->getNumberOfNodes() - 1;
      while(!stack.empty()) {
        int nIdx = stack.top();
        stack.pop();
        preorder2[count] = nIdx;
        count--;
        depth2 = std::max((int)predecessors2[nIdx].size(), depth2);
        std::vector<ftm::idNode> children;
        tree2->getChildren(nIdx, children);
        for(int cIdx : children) {
          stack.push(cIdx);
          predecessors2[cIdx].reserve(predecessors2[nIdx].size() + 1);
          predecessors2[cIdx].insert(predecessors2[cIdx].end(),
                                     predecessors2[nIdx].begin(),
                                     predecessors2[nIdx].end());
          predecessors2[cIdx].push_back(nIdx);
        }
      }

      // initialize memoization tables

      size_t nn1 = tree1->getNumberOfNodes();
      size_t nn2 = tree2->getNumberOfNodes();
      size_t dim1 = 1;
      size_t dim2 = (nn1 + 1) * dim1;
      size_t dim3 = (depth1 + 1) * dim2;
      size_t dim4 = (nn2 + 1) * dim3;

      //std::cout << (nn1 + 1) * (depth1 + 1) * (nn2 + 1) * (depth2 + 1) * sizeof(dataType) << std::endl;
      std::vector<dataType> memT((nn1 + 1) * (depth1 + 1) * (nn2 + 1)
                                 * (depth2 + 1));

      memT[nn1 + 0 * dim2 + nn2 * dim3 + 0 * dim4] = 0;
      for(size_t i = 0; i < nn1; i++) {
        int curr1 = preorder1[i];
        std::vector<ftm::idNode> children1;
        tree1->getChildren(curr1, children1);
        for(size_t l = 1; l <= predecessors1[preorder1[i]].size(); l++) {
          int parent1 = predecessors1[preorder1[i]]
                                     [predecessors1[preorder1[i]].size() - l];

          //-----------------------------------------------------------------------
          // Delete curr path and full subtree rooted in path
          memT[curr1 + l * dim2 + nn2 * dim3 + 0 * dim4]
            = editCost_Persistence<dataType>(
              curr1, parent1, -1, -1, tree1, tree2);
          for(auto child1 : children1) {
            memT[curr1 + l * dim2 + nn2 * dim3 + 0 * dim4]
              += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
          }
        }
      }
      for(size_t j = 0; j < nn2; j++) {
        int curr2 = preorder2[j];
        std::vector<ftm::idNode> children2;
        tree2->getChildren(curr2, children2);
        for(size_t l = 1; l <= predecessors2[preorder2[j]].size(); l++) {
          int parent2 = predecessors2[preorder2[j]]
                                     [predecessors2[preorder2[j]].size() - l];

          //-----------------------------------------------------------------------
          // Delete curr path and full subtree rooted in path
          memT[nn1 + 0 * dim2 + curr2 * dim3 + l * dim4]
            = editCost_Persistence<dataType>(
              -1, -1, curr2, parent2, tree1, tree2);
          for(auto child2 : children2) {
            memT[nn1 + 0 * dim2 + curr2 * dim3 + l * dim4]
              += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
          }
        }
      }

      for(size_t i = 0; i < nn1; i++) {
        int curr1 = preorder1[i];
        std::vector<ftm::idNode> children1;
        tree1->getChildren(curr1, children1);
        for(size_t j = 0; j < nn2; j++) {
          int curr2 = preorder2[j];
          std::vector<ftm::idNode> children2;
          tree2->getChildren(curr2, children2);
          for(size_t l1 = 1; l1 <= predecessors1[preorder1[i]].size(); l1++) {
            int parent1
              = predecessors1[preorder1[i]]
                             [predecessors1[preorder1[i]].size() - l1];
            for(size_t l2 = 1; l2 <= predecessors2[preorder2[j]].size(); l2++) {
              int parent2
                = predecessors2[preorder2[j]]
                               [predecessors2[preorder2[j]].size() - l2];

              //===============================================================================
              // If both trees not empty, find optimal edit operation

              //---------------------------------------------------------------------------
              // If both trees only have one branch, return edit cost between
              // the two branches
              if(tree1->getNumberOfChildren(curr1) == 0
                 and tree2->getNumberOfChildren(curr2) == 0) {
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4]
                  = editCost_Persistence<dataType>(
                    curr1, parent1, curr2, parent2, tree1, tree2);
              }
              //---------------------------------------------------------------------------
              // If first tree only has one branch, try all decompositions of
              // second tree
              else if(tree1->getNumberOfChildren(curr1) == 0) {
                dataType d = std::numeric_limits<dataType>::max();
                for(auto child2_mb : children2) {
                  dataType d_ = memT[curr1 + l1 * dim2 + child2_mb * dim3
                                     + (l2 + 1) * dim4];
                  for(auto child2 : children2) {
                    if(child2 == child2_mb) {
                      continue;
                    }
                    d_ += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
                  }
                  d = std::min(d, d_);
                }
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] = d;
              }
              //---------------------------------------------------------------------------
              // If second tree only has one branch, try all decompositions of
              // first tree
              else if(tree2->getNumberOfChildren(curr2) == 0) {
                dataType d = std::numeric_limits<dataType>::max();
                for(auto child1_mb : children1) {
                  dataType d_ = memT[child1_mb + (l1 + 1) * dim2 + curr2 * dim3
                                     + l2 * dim4];
                  for(auto child1 : children1) {
                    if(child1 == child1_mb) {
                      continue;
                    }
                    d_ += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
                  }
                  d = std::min(d, d_);
                }
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] = d;
              }
              //---------------------------------------------------------------------------
              // If both trees have more than one branch, try all decompositions
              // of both trees
              else {
                dataType d = std::numeric_limits<dataType>::max();
                //-----------------------------------------------------------------------
                // Try all possible main branches of first tree (child1_mb) and
                // all possible main branches of second tree (child2_mb) Then
                // try all possible matchings of subtrees
                if(tree1->getNumberOfChildren(curr1) == 2
                   && tree2->getNumberOfChildren(curr2) == 2) {
                  int child11 = children1[0];
                  int child12 = children1[1];
                  int child21 = children2[0];
                  int child22 = children2[1];
                  d = std::min<dataType>(
                    d, memT[child11 + 1 * dim2 + child21 * dim3 + 1 * dim4]
                         + memT[child12 + 1 * dim2 + child22 * dim3 + 1 * dim4]
                         + editCost_Persistence<dataType>(
                           curr1, parent1, curr2, parent2, tree1, tree2));
                  d = std::min<dataType>(
                    d, memT[child11 + 1 * dim2 + child22 * dim3 + 1 * dim4]
                         + memT[child12 + 1 * dim2 + child21 * dim3 + 1 * dim4]
                         + editCost_Persistence<dataType>(
                           curr1, parent1, curr2, parent2, tree1, tree2));
                } else {
                  auto f = [&](int r, int c) {
                    size_t c1 = r < tree1->getNumberOfChildren(curr1)
                                  ? children1[r]
                                  : nn1;
                    size_t c2 = c < tree2->getNumberOfChildren(curr2)
                                  ? children2[c]
                                  : nn2;
                    int l1_ = c1 == nn1 ? 0 : 1;
                    int l2_ = c2 == nn2 ? 0 : 1;
                    return memT[c1 + l1_ * dim2 + c2 * dim3 + l2_ * dim4];
                  };
                  int size = std::max(tree1->getNumberOfChildren(curr1),
                                      tree2->getNumberOfChildren(curr2))
                             + 1;
                  auto costMatrix = std::vector<std::vector<dataType>>(
                    size, std::vector<dataType>(size, 0));
                  std::vector<MatchingType> matching;
                  for(int r = 0; r < size; r++) {
                    for(int c = 0; c < size; c++) {
                      costMatrix[r][c] = f(r, c);
                    }
                  }

                  AssignmentSolver<dataType> *assignmentSolver;
                  AssignmentExhaustive<dataType> solverExhaustive;
                  AssignmentMunkres<dataType> solverMunkres;
                  AssignmentAuction<dataType> solverAuction;
                  switch(assignmentSolverID_) {
                    case 1:
                      solverExhaustive = AssignmentExhaustive<dataType>();
                      assignmentSolver = &solverExhaustive;
                      break;
                    case 2:
                      solverMunkres = AssignmentMunkres<dataType>();
                      assignmentSolver = &solverMunkres;
                      break;
                    case 0:
                    default:
                      solverAuction = AssignmentAuction<dataType>();
                      assignmentSolver = &solverAuction;
                  }
                  assignmentSolver->setInput(costMatrix);
                  assignmentSolver->setBalanced(true);
                  assignmentSolver->run(matching);
                  dataType d_ = editCost_Persistence<dataType>(
                    curr1, parent1, curr2, parent2, tree1, tree2);
                  for(auto m : matching)
                    d_ += std::get<2>(m);
                  d = std::min(d, d_);
                }
                //-----------------------------------------------------------------------
                // Try to continue main branch on one child of first tree and
                // delete all other subtrees Then match continued branch to
                // current branch in second tree
                for(auto child1_mb : children1) {
                  dataType d_ = memT[child1_mb + (l1 + 1) * dim2 + curr2 * dim3
                                     + l2 * dim4];
                  for(auto child1 : children1) {
                    if(child1 == child1_mb) {
                      continue;
                    }
                    d_ += memT[child1 + 1 * dim2 + nn2 * dim3 + 0 * dim4];
                  }
                  d = std::min(d, d_);
                }
                //-----------------------------------------------------------------------
                // Try to continue main branch on one child of second tree and
                // delete all other subtrees Then match continued branch to
                // current branch in first tree
                for(auto child2_mb : children2) {
                  dataType d_ = memT[curr1 + l1 * dim2 + child2_mb * dim3
                                     + (l2 + 1) * dim4];
                  for(auto child2 : children2) {
                    if(child2 == child2_mb) {
                      continue;
                    }
                    d_ += memT[nn1 + 0 * dim2 + child2 * dim3 + 1 * dim4];
                  }
                  d = std::min(d, d_);
                }
                memT[curr1 + l1 * dim2 + curr2 * dim3 + l2 * dim4] = d;
              }
            }
          }
        }
      }

      std::vector<ftm::idNode> children1;
      tree1->getChildren(rootID1, children1);
      std::vector<ftm::idNode> children2;
      tree2->getChildren(rootID2, children2);

      dataType res
        = memT[children1[0] + 1 * dim2 + children2[0] * dim3 + 1 * dim4];

      if(computeMapping_ && outputMatching){

        outputMatching->clear();
        traceMapping_path(tree1,tree2,children1[0],1,children2[0],1,predecessors1,predecessors2,depth1,depth2,memT,*outputMatching);

        // dataType cost_mapping = 0;
        // dataType cost_ins = 0;
        // dataType cost_del = 0;
        // std::vector<bool> matchedNodes1(tree1->getNumberOfNodes(),false);
        // std::vector<bool> matchedNodes2(tree2->getNumberOfNodes(),false);
        // for(auto m : *outputMatching){
        //   matchedNodes1[m.first.first] = true;
        //   matchedNodes1[m.first.second] = true;
        //   matchedNodes2[m.second.first] = true;
        //   matchedNodes2[m.second.second] = true;
        //   ftm::idNode cn = tree1->getParentSafe(m.first.first);
        //   while(cn!=m.first.second){
        //     matchedNodes1[cn] = true;
        //     cn = tree1->getParentSafe(cn);
        //   }
        //   cn = tree2->getParentSafe(m.second.first);
        //   while(cn!=m.second.second){
        //     matchedNodes2[cn] = true;
        //     cn = tree2->getParentSafe(cn);
        //   }
        //   dataType cost = editCost_Persistence<dataType>(
        //                     m.first.first, m.first.second, m.second.first, m.second.second, tree1, tree2);
        //   cost_mapping += cost;
        //   // std::cout << "   (" << m.first.first << " " << m.first.second << ") - (" << m.second.first << " " << m.second.second << ") : " << cost << "\n";
        // }
        // for(ftm::idNode i=0; i<tree1->getNumberOfNodes(); i++){
        //   if(!matchedNodes1[i]){
        //     dataType cost = editCost_Persistence<dataType>(i, tree1->getParentSafe(i), -1, -1, tree1, tree2);
        //     // std::cout << "   (" << i << " " << tree1->getParentSafe(i) << ") - (" << -1 << " " << -1 << ") : " << cost << "\n";
        //     cost_del += cost;
        //   }
        // }
        // for(ftm::idNode i=0; i<tree2->getNumberOfNodes(); i++){
        //   if(!matchedNodes2[i]){
        //     dataType cost = editCost_Persistence<dataType>(-1, -1, i, tree2->getParentSafe(i), tree1, tree2);
        //     // std::cout << "   (" << -1 << " " << -1 << ") - (" << i << " " << tree2->getParentSafe(i) << ") : " << cost << "\n";
        //     cost_ins += cost;
        //   }
        // }
        // // std::cout << res << " " << cost_mapping+cost_ins+cost_del << std::endl;
        // // std::cout << res << " " << cost_mapping << std::endl;

      }

      return squared_ ? std::sqrt(res) : res;
    }

    template <class dataType>
    dataType editDistance_path(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2, std::vector<std::tuple<ftm::idNode, ftm::idNode, double>> *outputMatching){
      
      std::vector<int> matchedNodes(tree1->getNumberOfNodes(),-1);
      std::vector<std::pair<std::pair<ftm::idNode, ftm::idNode>,std::pair<ftm::idNode, ftm::idNode>>> mapping;
      dataType res = editDistance_path<dataType>(tree1,tree2,&mapping);
      if(computeMapping_ && outputMatching){
        outputMatching->clear();
        for(auto m : mapping){
          matchedNodes[m.first.first] = m.second.first;
          matchedNodes[m.first.second] = m.second.second;
        }
        for(ftm::idNode i=0; i<matchedNodes.size(); i++){
          if(matchedNodes[i]>=0) outputMatching->emplace_back(std::make_tuple(i,matchedNodes[i], 0.0));
        }
      }

      return res;

    }

    template <class dataType>
    dataType editDistance_path(ftm::FTMTree_MT *tree1, ftm::FTMTree_MT *tree2){
      return editDistance_path<dataType>(tree1,tree2,(std::vector<std::pair<std::pair<ftm::idNode, ftm::idNode>,std::pair<ftm::idNode, ftm::idNode>>>*) nullptr);
    }
  };
} // namespace ttk
