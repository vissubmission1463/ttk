/// \ingroup vtk
/// \class ttkMergeTreeVisualization
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Visualization module for merge trees.

#pragma once

#include <FTMTree.h>
#include <MergeTreeVisualization.h>

#include <ttkAlgorithm.h>

// VTK Includes
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

class ttkMergeTreeVisualization : public ttk::MergeTreeVisualization {
private:
  // redefine types
  using idNode = ttk::ftm::idNode;
  using SimplexId = ttk::SimplexId;
  using FTMTree_MT = ttk::ftm::FTMTree_MT;

  // Visualization parameters
  bool PlanarLayout = false;
  double DimensionSpacing = 1.;
  int DimensionToShift = 0;
  double XShift = 1.0;
  double YShift = 0.0;
  double ZShift = 0.0;
  bool OutputSegmentation = false;
  int MaximumImportantPairs = 0;
  int MinimumImportantPairs = 0;
  bool outputTreeNodeIndex = false;
  bool isPersistenceDiagram = false;
  bool isPDSadMax = true;
  bool enableBarycenterAlignment = false;

  // Shift mode
  // -1: None ; 0: Star ; 1: Star Barycenter ; 2: Line ; 3: Double Line
  int ShiftMode = 0;

  // Offset
  int iSampleOffset = 0;
  int noSampleOffset = 0;
  double prevXMaxOffset = 0;

  // Print only one tree
  int printTreeId = -1; // -1 for all
  int printClusterId = -1; // -1 for all

  // Barycenter position according alpha
  bool BarycenterPositionAlpha = false;
  double Alpha = 0.5;

  // Used for critical type and for point coordinates
  std::vector<vtkUnstructuredGrid *> treesNodes;
  std::vector<std::vector<int>>
    treesNodeCorrMesh; // used to acces treesNodes given input trees

  // Segmentation
  std::vector<vtkDataSet *> treesSegmentation;

  // Clustering output
  std::vector<int> clusteringAssignment;
  std::vector<std::vector<std::vector<std::tuple<idNode, idNode, double>>>>
    outputMatchingBarycenter;

  // Path layout
  std::vector<std::vector<
    std::vector<std::pair<std::pair<ttk::ftm::idNode, ttk::ftm::idNode>,
                          std::pair<ttk::ftm::idNode, ttk::ftm::idNode>>>>>
    pathMatchings;

  // Barycenter output
  std::vector<std::vector<float>> allBaryPercentMatch;

  // Temporal Subsampling Output
  std::vector<bool> interpolatedTrees;

  // Output
  vtkUnstructuredGrid *vtkOutputNode{};
  vtkUnstructuredGrid *vtkOutputArc{};
  vtkDataSet *vtkOutputSegmentation{};

  // Matching output
  vtkUnstructuredGrid *vtkOutputNode1{}, *vtkOutputNode2{}; // input data
  std::vector<std::vector<SimplexId>> nodeCorr1, nodeCorr2;
  vtkUnstructuredGrid *vtkOutputMatching{}; // output

  // Custom array
  std::vector<std::tuple<std::string, std::vector<double>>> customArrays;
  std::vector<std::tuple<std::string, std::vector<int>>> customIntArrays;
  std::vector<std::tuple<std::string, std::vector<std::string>>>
    customStringArrays;

  // Filled by the algorithm
  std::vector<std::vector<SimplexId>> nodeCorr;
  std::vector<double> clusterShift;
  double prevXMax = 0;

public:
  ttkMergeTreeVisualization() = default;
  ;
  ~ttkMergeTreeVisualization() override = default;
  ;

  // ==========================================================================
  // Getter / Setter
  // ==========================================================================
  // Visualization parameters
  void setPlanarLayout(bool b) {
    PlanarLayout = b;
  }
  void setDimensionSpacing(double d) {
    DimensionSpacing = d;
  }
  void setDimensionToShift(int i) {
    DimensionToShift = i;
  }
  void setXshift(double shift) {
    XShift = shift;
  }
  void setYshift(double shift) {
    YShift = shift;
  }
  void setZshift(double shift) {
    ZShift = shift;
  }
  void setDimensionsShift(double xShift, double yShift, double zShift) {
    setXshift(xShift);
    setYshift(yShift);
    setZshift(zShift);
  }
  void setOutputSegmentation(bool b) {
    OutputSegmentation = b;
  }
  void setShiftMode(int mode) {
    ShiftMode = mode;
  }
  void setMaximumImportantPairs(int maxPairs) {
    MaximumImportantPairs = maxPairs;
  }
  void setMinimumImportantPairs(int minPairs) {
    MinimumImportantPairs = minPairs;
  }

  void setOutputTreeNodeId(int doOutput) {
    outputTreeNodeIndex = doOutput;
  }

  void setIsPersistenceDiagram(bool isPD) {
    isPersistenceDiagram = isPD;
  }
  void setIsPDSadMax(bool isSadMax) {
    isPDSadMax = isSadMax;
  }

  void setEnableBarycenterAlignment(bool eba) {
    enableBarycenterAlignment = eba;
  }

  // Offset
  void setISampleOffset(int offset) {
    iSampleOffset = offset;
  }
  void setNoSampleOffset(int offset) {
    noSampleOffset = offset;
  }
  void setPrevXMaxOffset(double offset) {
    prevXMaxOffset = offset;
  }

  // Print only one tree
  void setPrintTreeId(int id) {
    printTreeId = id;
  }
  void setPrintClusterId(int id) {
    printClusterId = id;
  }

  // Barycenter position according alpha
  void setBarycenterPositionAlpha(bool pos) {
    BarycenterPositionAlpha = pos;
  }
  void setAlpha(double alpha) {
    Alpha = alpha;
  }

  // Used for critical type and if not planar layout for point coordinates
  void setTreesNodes(std::vector<vtkUnstructuredGrid *> &nodes) {
    treesNodes = nodes;
  }
  void setTreesNodeCorrMesh(std::vector<std::vector<int>> &nodeCorrMesh) {
    treesNodeCorrMesh = nodeCorrMesh;
  }
  void setTreesNodes(vtkUnstructuredGrid *nodes) {
    treesNodes.clear();
    treesNodes.emplace_back(nodes);
  }
  void setTreesNodeCorrMesh(std::vector<int> &nodeCorrMesh) {
    treesNodeCorrMesh.clear();
    treesNodeCorrMesh.emplace_back(nodeCorrMesh);
  }

  // Segmentation
  void setTreesSegmentation(std::vector<vtkDataSet *> &segmentation) {
    treesSegmentation = segmentation;
  }
  void setTreesSegmentation(vtkDataSet *segmentation) {
    treesSegmentation.clear();
    treesSegmentation.emplace_back(segmentation);
  }

  // Clustering output
  void setClusteringAssignment(std::vector<int> &asgn) {
    clusteringAssignment = asgn;
  }
  void setOutputMatchingBarycenter(
    std::vector<std::vector<std::vector<std::tuple<idNode, idNode, double>>>>
      &matching) {
    outputMatchingBarycenter = matching;
  }

  // Path layout
  void setPathMatchings(
    std::vector<std::vector<
      std::vector<std::pair<std::pair<ttk::ftm::idNode, ttk::ftm::idNode>,
                            std::pair<ttk::ftm::idNode, ttk::ftm::idNode>>>>>
      &matchings) {
    pathMatchings = matchings;
  }
  void getTreePathing(
    FTMTree_MT *tree,
    std::vector<
      std::vector<std::pair<std::pair<ttk::ftm::idNode, ttk::ftm::idNode>,
                            std::pair<ttk::ftm::idNode, ttk::ftm::idNode>>>>
      &matchings,
    bool isFirstTree,
    std::vector<ttk::ftm::idNode> &pathing,
    std::vector<int> &pathingID) {
    pathing = std::vector<idNode>(tree->getNumberOfNodes());
    pathingID = std::vector<int>(tree->getNumberOfNodes(), -1);
    std::vector<int> nodeLevel;
    tree->getAllNodeLevel(nodeLevel);
    int pathID = 0;
    std::vector<std::vector<bool>> processed(
      tree->getNumberOfNodes(),
      std::vector<bool>(tree->getNumberOfNodes(), false));
    for(auto &matching : matchings) {
      for(auto &match : matching) {
        auto first = (isFirstTree ? match.first.first : match.second.first);
        auto second = (isFirstTree ? match.first.second : match.second.second);
        auto lowest = (nodeLevel[first] < nodeLevel[second] ? second : first);
        auto highest = (nodeLevel[first] < nodeLevel[second] ? first : second);
        if(processed[lowest][highest])
          continue;
        processed[lowest][highest] = true;
        while(lowest != highest) {
          pathing[lowest] = highest;
          pathingID[lowest] = pathID;
          lowest = tree->getParentSafe(lowest);
        }
        if(tree->isRoot(highest)) {
          pathing[highest] = highest;
          pathingID[highest] = pathID;
        }
        ++pathID;
      }
    }
  }
  void getTreePathing(
    FTMTree_MT *tree,
    std::vector<std::pair<std::pair<ttk::ftm::idNode, ttk::ftm::idNode>,
                          std::pair<ttk::ftm::idNode, ttk::ftm::idNode>>>
      &matching,
    bool isFirstTree,
    std::vector<ttk::ftm::idNode> &pathing,
    std::vector<int> &pathingID) {
    std::vector<
      std::vector<std::pair<std::pair<ttk::ftm::idNode, ttk::ftm::idNode>,
                            std::pair<ttk::ftm::idNode, ttk::ftm::idNode>>>>
      matchings{matching};
    getTreePathing(tree, matchings, isFirstTree, pathing, pathingID);
  }

  // Barycenter output
  std::vector<std::vector<float>> getAllBaryPercentMatch() {
    return allBaryPercentMatch;
  }
  void
    setAllBaryPercentMatch(std::vector<std::vector<float>> &baryPercentMatch) {
    allBaryPercentMatch = baryPercentMatch;
  }

  // Temporal Subsampling Output
  void setInterpolatedTrees(std::vector<bool> &isInterpolatedTrees) {
    interpolatedTrees = isInterpolatedTrees;
  }

  // Output
  void setVtkOutputNode(vtkUnstructuredGrid *vtkNode) {
    vtkOutputNode = vtkNode;
  }
  void setVtkOutputArc(vtkUnstructuredGrid *vtkArc) {
    vtkOutputArc = vtkArc;
  }
  void setVtkOutputSegmentation(vtkDataSet *vtkSegmentation) {
    vtkOutputSegmentation = vtkSegmentation;
  }

  // Matching output
  void setVtkOutputNode1(vtkUnstructuredGrid *vtkNode1) {
    vtkOutputNode1 = vtkNode1;
  }
  void setVtkOutputNode2(vtkUnstructuredGrid *vtkNode2) {
    vtkOutputNode2 = vtkNode2;
  }
  void setNodeCorr1(std::vector<std::vector<SimplexId>> &nodeCorrT) {
    nodeCorr1 = nodeCorrT;
  }
  void setNodeCorr2(std::vector<std::vector<SimplexId>> &nodeCorrT) {
    nodeCorr2 = nodeCorrT;
  }
  void setVtkOutputMatching(vtkUnstructuredGrid *vtkMatching) {
    vtkOutputMatching = vtkMatching;
  }
  void setOutputMatching(
    std::vector<std::tuple<idNode, idNode, double>> &matching) {
    outputMatchingBarycenter.resize(1);
    outputMatchingBarycenter[0].resize(1);
    outputMatchingBarycenter[0][0] = matching;
  }

  // Custom array
  void addCustomArray(std::string &name, std::vector<double> &vec) {
    customArrays.emplace_back(name, vec);
  }
  void addCustomIntArray(std::string &name, std::vector<int> &vec) {
    customIntArrays.emplace_back(name, vec);
  }
  void addCustomStringArray(std::string &name, std::vector<std::string> &vec) {
    customStringArrays.emplace_back(name, vec);
  }
  void clearCustomArrays() {
    customArrays.clear();
  }
  void clearCustomIntArrays() {
    customIntArrays.clear();
  }
  void clearCustomStringArrays() {
    customStringArrays.clear();
  }
  void clearAllCustomArrays() {
    clearCustomArrays();
    clearCustomIntArrays();
    clearCustomStringArrays();
  }

  template <class dataType>
  void addVtkCustomArrays(
    std::vector<std::tuple<std::string, std::vector<dataType>>> &cArrays,
    std::vector<std::vector<dataType>> &cArraysValues,
    vtkUnstructuredGrid *vtkOutput,
    int type,
    int output) {
    for(unsigned int i = 0; i < cArrays.size(); ++i) {
      vtkNew<vtkDoubleArray> customDoubleArrayVtk;
      vtkNew<vtkIntArray> customIntArrayVtk;
      vtkNew<vtkStringArray> customStringArrayVtk;
      vtkAbstractArray *customArrayVtk;
      if(type == 0)
        customArrayVtk = customDoubleArrayVtk;
      else if(type == 1)
        customArrayVtk = customIntArrayVtk;
      else
        customArrayVtk = customStringArrayVtk;
      customArrayVtk->SetName(std::get<0>(cArrays[i]).c_str());
      customArrayVtk->SetNumberOfTuples(cArraysValues[i].size());
      for(unsigned int j = 0; j < cArraysValues[i].size(); ++j) {
        // Add value depending on type (vtkAbstractArray can not be used here)
        if(type == 0) {
          double doubleValue = (*(std::vector<double> *)&(cArraysValues[i]))[j];
          customDoubleArrayVtk->SetValue(j, doubleValue);
        } else if(type == 1) {
          int intValue = (*(std::vector<int> *)&(cArraysValues[i]))[j];
          customIntArrayVtk->SetValue(j, intValue);
        } else {
          std::string stringValue
            = (*(std::vector<std::string> *)&(cArraysValues[i]))[j];
          customStringArrayVtk->SetValue(j, stringValue);
        }
      }
      if(output == 0)
        vtkOutput->GetPointData()->AddArray(customArrayVtk);
      else
        vtkOutput->GetCellData()->AddArray(customArrayVtk);
    }
  }

  void getTreeNodeIdRev(vtkDataArray *treeNodeIdArray,
                        std::vector<int> &treeNodeIdRev) {
    double valueRange[2];
    treeNodeIdArray->GetRange(valueRange);
    int maxValue = valueRange[1];
    treeNodeIdRev.clear();
    treeNodeIdRev.resize(maxValue + 1);
    for(int i = 0; i < treeNodeIdArray->GetNumberOfValues(); ++i)
      treeNodeIdRev[treeNodeIdArray->GetTuple1(i)] = i;
  }

  void copyPointData(vtkUnstructuredGrid *treeNodes,
                     std::vector<int> &nodeCorrT) {
    if(!treeNodes)
      return;

    auto treeNodeIdArray = treeNodes->GetPointData()->GetArray("TreeNodeId");
    std::vector<int> treeNodeIdRev;
    if(treeNodeIdArray)
      getTreeNodeIdRev(treeNodeIdArray, treeNodeIdRev);

    for(int i = 0; i < treeNodes->GetPointData()->GetNumberOfArrays(); ++i) {
      auto dataArray
        = vtkDataArray::SafeDownCast(treeNodes->GetPointData()->GetArray(i));
      auto stringArray
        = vtkStringArray::SafeDownCast(treeNodes->GetPointData()->GetArray(i));
      vtkAbstractArray *array;
      if(dataArray)
        array = dataArray;
      else if(stringArray)
        array = stringArray;
      else
        continue;
      auto vecSize = (nodeCorrT.size() == 0 ? array->GetNumberOfValues()
                                            : nodeCorrT.size());
      std::vector<double> vec(vecSize);
      std::vector<std::string> vecString(vecSize);
      for(unsigned int j = 0; j < vec.size(); ++j) {
        int toGet = (nodeCorrT.size() == 0 ? j : nodeCorrT[j]);
        if(treeNodeIdArray)
          toGet = (nodeCorrT.size() == 0 ? treeNodeIdRev[j]
                                         : treeNodeIdRev[nodeCorrT[j]]);
        auto value = array->GetVariantValue(toGet);
        if(dataArray)
          vec[j] = value.ToDouble();
        else
          vecString[j] = value.ToString();
      }
      std::string name{array->GetName()};
      if(dataArray)
        addCustomArray(name, vec);
      else
        addCustomStringArray(name, vecString);
    }
  }
  void copyPointData(vtkUnstructuredGrid *treeNodes) {
    std::vector<int> nodeCorrT;
    copyPointData(treeNodes, nodeCorrT);
  }

  // Filled by the algorithm
  std::vector<std::vector<SimplexId>> getNodeCorr() {
    return nodeCorr;
  }
  std::vector<double> getClusterShift() {
    return clusterShift;
  }
  double getPrevXMax() {
    return prevXMax;
  }

  // ==========================================================================
  // Matching Visualization
  // ==========================================================================
  template <class dataType>
  void makeMatchingOutput(FTMTree_MT *tree1, FTMTree_MT *tree2) {
    std::vector<FTMTree_MT *> trees{tree1, tree2};
    std::vector<FTMTree_MT *> barycenters;

    makeMatchingOutput<dataType>(trees, barycenters);
  }

  template <class dataType>
  void makeMatchingOutput(std::vector<FTMTree_MT *> &trees,
                          std::vector<FTMTree_MT *> &barycenters) {
    int numInputs = trees.size();
    int NumberOfBarycenters = barycenters.size();
    bool clusteringOutput = (NumberOfBarycenters != 0);
    NumberOfBarycenters
      = std::max(NumberOfBarycenters, 1); // to always enter the outer loop
    if(not clusteringOutput)
      numInputs = 1;

    vtkNew<vtkUnstructuredGrid> vtkMatching{};
    vtkNew<vtkPoints> pointsM{};

    // Fields
    vtkNew<vtkIntArray> matchingID{};
    matchingID->SetName("MatchingID");
    vtkNew<vtkIntArray> matchingType{};
    matchingType->SetName("MatchingType");
    vtkNew<vtkDoubleArray> matchPers{};
    matchPers->SetName("MeanMatchedPersistence");
    vtkNew<vtkDoubleArray> costArray{};
    costArray->SetName("Cost");
    vtkNew<vtkIntArray> tree1NodeIdField{};
    tree1NodeIdField->SetName("tree1NodeId");
    vtkNew<vtkIntArray> tree2NodeIdField{};
    tree2NodeIdField->SetName("tree2NodeId");
    vtkNew<vtkIntArray> mergeTree1NodeIdField{};
    mergeTree1NodeIdField->SetName("mergeTree1NodeId");
    vtkNew<vtkIntArray> mergeTree2NodeIdField{};
    mergeTree2NodeIdField->SetName("mergeTree2NodeId");
    vtkNew<vtkIntArray> isBarycenterNodeField{};
    isBarycenterNodeField->SetName("isBarycenterNode");

    vtkNew<vtkFloatArray> matchingPercentMatch{};
    matchingPercentMatch->SetName("MatchingPercentMatch");

    // Iterate through clusters and trees
    printMsg(
      "// Iterate through clusters and trees", ttk::debug::Priority::VERBOSE);
    int count = 0;
    for(int c = 0; c < NumberOfBarycenters; ++c) {
      for(int i = 0; i < numInputs; ++i) {
        if((printTreeId == -1 and printClusterId != -1 and c != printClusterId)
           or (printTreeId != -1 and printClusterId == -1 and i != printTreeId)
           or (printTreeId != -1 and printClusterId != -1
               and (c != printClusterId or i != printTreeId)))
          continue;
        for(std::tuple<idNode, idNode, double> match :
            outputMatchingBarycenter[c][i]) {
          vtkIdType pointIds[2];
          idNode tree1NodeId = std::get<0>(match);
          idNode tree2NodeId = std::get<1>(match);
          double cost = std::get<2>(match);
          FTMTree_MT *tree1;
          FTMTree_MT *tree2 = trees[i];
          if(not clusteringOutput) {
            tree1 = trees[0];
            tree2 = trees[1];
          } else
            tree1 = barycenters[c];

          // Get first point
          printMsg("// Get first point", ttk::debug::Priority::VERBOSE);
          SimplexId pointToGet1 = clusteringOutput ? nodeCorr2[c][tree1NodeId]
                                                   : nodeCorr1[0][tree1NodeId];
          double *point1 = vtkOutputNode2->GetPoints()->GetPoint(pointToGet1);
          const SimplexId nextPointId1 = pointsM->InsertNextPoint(point1);
          if(not clusteringOutput) isBarycenterNodeField->InsertNextTuple1(0);
          else isBarycenterNodeField->InsertNextTuple1(1);
          pointIds[0] = nextPointId1;

          // Get second point
          printMsg("// Get second point", ttk::debug::Priority::VERBOSE);
          SimplexId pointToGet2 = clusteringOutput ? nodeCorr1[i][tree2NodeId]
                                                   : nodeCorr1[1][tree2NodeId];
          double *point2 = vtkOutputNode1->GetPoints()->GetPoint(pointToGet2);
          const SimplexId nextPointId2 = pointsM->InsertNextPoint(point2);
          isBarycenterNodeField->InsertNextTuple1(0);
          pointIds[1] = nextPointId2;

          // Add cell
          printMsg("// Add cell", ttk::debug::Priority::VERBOSE);
          vtkMatching->InsertNextCell(VTK_LINE, 2, pointIds);

          // Add arc matching percentage
          printMsg(
            "// Add arc matching percentage", ttk::debug::Priority::VERBOSE);
          if(allBaryPercentMatch.size() != 0)
            matchingPercentMatch->InsertNextTuple1(
              allBaryPercentMatch[c][tree1NodeId]);

          // Add tree1 and tree2 node ids
          printMsg(
            "// Add tree1 and tree2 node ids", ttk::debug::Priority::VERBOSE);
          tree1NodeIdField->InsertNextTuple1(pointToGet1);
          tree2NodeIdField->InsertNextTuple1(pointToGet2);
          mergeTree1NodeIdField->InsertNextTuple1(tree1NodeId);
          mergeTree2NodeIdField->InsertNextTuple1(tree2NodeId);

          // Add matching ID
          matchingID->InsertNextTuple1(count);

          // Add matching type
          printMsg("// Add matching type", ttk::debug::Priority::VERBOSE);
          int thisType = 0;
          int tree1NodeDown
            = tree1->getNode(tree1NodeId)->getNumberOfDownSuperArcs();
          int tree1NodeUp
            = tree1->getNode(tree1NodeId)->getNumberOfUpSuperArcs();
          int tree2NodeDown
            = tree2->getNode(tree2NodeId)->getNumberOfDownSuperArcs();
          int tree2NodeUp
            = tree2->getNode(tree2NodeId)->getNumberOfUpSuperArcs();
          if(tree1NodeDown != 0 and tree1NodeUp != 0 and tree2NodeDown != 0
             and tree2NodeUp != 0)
            thisType = 1; // Saddle to Saddle
          if(tree1NodeDown == 0 and tree1NodeUp != 0 and tree2NodeDown == 0
             and tree2NodeUp != 0)
            thisType = 2; // Leaf to leaf
          if(tree1NodeDown != 0 and tree1NodeUp == 0 and tree2NodeDown != 0
             and tree2NodeUp == 0)
            thisType = 3; // Root to root
          matchingType->InsertNextTuple1(thisType);

          // Add mean matched persistence
          printMsg(
            "// Add mean matched persistence", ttk::debug::Priority::VERBOSE);
          double tree1Pers = tree1->getNodePersistence<dataType>(tree1NodeId);
          double tree2Pers = tree2->getNodePersistence<dataType>(tree2NodeId);
          double meanPersistence = (tree1Pers + tree2Pers) / 2;
          matchPers->InsertNextTuple1(meanPersistence);

          // Add cost
          costArray->InsertNextTuple1(cost);

          count++;
        }
      }
    }
    vtkMatching->SetPoints(pointsM);
    vtkMatching->GetCellData()->AddArray(matchingType);
    vtkMatching->GetCellData()->AddArray(matchPers);
    vtkMatching->GetCellData()->AddArray(matchingID);
    vtkMatching->GetCellData()->AddArray(costArray);
    vtkMatching->GetCellData()->AddArray(tree1NodeIdField);
    vtkMatching->GetCellData()->AddArray(tree2NodeIdField);
    vtkMatching->GetCellData()->AddArray(mergeTree1NodeIdField);
    vtkMatching->GetCellData()->AddArray(mergeTree2NodeIdField);
    vtkMatching->GetPointData()->AddArray(isBarycenterNodeField);
    if(allBaryPercentMatch.size() != 0)
      vtkMatching->GetCellData()->AddArray(matchingPercentMatch);
    vtkOutputMatching->ShallowCopy(vtkMatching);
  }

  // ==========================================================================
  // Trees Visualization
  // ==========================================================================
  template <class dataType>
  void makeTreesOutput(FTMTree_MT *tree1) {
    std::vector<FTMTree_MT *> trees{tree1};

    makeTreesOutput<dataType>(trees);
  }

  template <class dataType>
  void makeTreesOutput(FTMTree_MT *tree1, FTMTree_MT *tree2) {
    std::vector<FTMTree_MT *> trees{tree1, tree2};

    makeTreesOutput<dataType>(trees);
  }

  template <class dataType>
  void makeTreesOutput(std::vector<FTMTree_MT *> &trees) {
    std::vector<FTMTree_MT *> barycenters;
    clusteringAssignment.clear();
    clusteringAssignment.resize(trees.size(), 0);

    makeTreesOutput<dataType>(trees, barycenters);
  }

  template <class dataType>
  void makeTreesOutput(std::vector<FTMTree_MT *> &trees,
                       std::vector<FTMTree_MT *> &barycenters) {
    int numInputs = trees.size();
    int numInputsOri = numInputs;
    int NumberOfBarycenters = barycenters.size();
    bool clusteringOutput = (NumberOfBarycenters != 0);
    NumberOfBarycenters
      = std::max(NumberOfBarycenters, 1); // to always enter the outer loop
    PlanarLayout |= isPersistenceDiagram;
    bool alignTrees = trees.size() == 2 and barycenters.size() == 1
                      and enableBarycenterAlignment;

    // TreeNodeIdRev
    for(int i = 0; i < numInputs; ++i) {
      if(i < (int)treesNodes.size() and treesNodes[i]) {
        auto treeNodeIdArray
          = treesNodes[i]->GetPointData()->GetArray("TreeNodeId");
        if(treeNodeIdArray) {
          std::vector<int> treeNodeIdRev;
          getTreeNodeIdRev(treeNodeIdArray, treeNodeIdRev);
          for(unsigned int j = 0; j < treesNodeCorrMesh[i].size(); ++j)
            treesNodeCorrMesh[i][j] = treeNodeIdRev[treesNodeCorrMesh[i][j]];
        }
      }
    }

    // Bounds
    printMsg("Bounds and branching", ttk::debug::Priority::VERBOSE);
    std::vector<std::tuple<double, double, double, double, double, double>>
      allBounds(numInputs);
    for(int i = 0; i < numInputs; ++i) {
      if(OutputSegmentation and treesSegmentation[i]) {
        double *tBounds = treesSegmentation[i]->GetBounds();
        allBounds[i] = std::make_tuple(tBounds[0], tBounds[1], tBounds[2],
                                       tBounds[3], tBounds[4], tBounds[5]);
      } else if(treesNodes.size() != 0 and treesNodes[i] != nullptr) {
        if(not isPersistenceDiagram)
          allBounds[i]
            = getRealBounds(treesNodes[i], trees[i], treesNodeCorrMesh[i]);
        else {
          double bounds[6];
          treesNodes[i]->GetBounds(bounds);
          allBounds[i] = std::make_tuple(
            bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
        }
      } else {
        allBounds[i] = allBounds[0];
      }
    }

    std::vector<std::tuple<double, double, double, double, double, double>>
      allBaryBounds(barycenters.size());
    std::vector<std::vector<idNode>> allBaryBranching(barycenters.size()),
      allBaryPathing(barycenters.size());
    std::vector<std::vector<int>> allBaryBranchingID(barycenters.size()),
      allBaryPathingID(barycenters.size());
    for(size_t c = 0; c < barycenters.size(); ++c) {
      allBaryBounds[c] = getMaximalBounds(allBounds, clusteringAssignment, c);
      if(not isPersistenceDiagram)
        barycenters[c]->getTreeBranching(
          allBaryBranching[c], allBaryBranchingID[c]);
      if(pathPlanarLayout_ and !pathMatchings.empty()
         and !pathMatchings[c].empty())
        getTreePathing(barycenters[c], pathMatchings[c], true,
                       allBaryPathing[c], allBaryPathingID[c]);
    }
    if(not clusteringOutput)
      allBaryBounds.emplace_back(
        getMaximalBounds(allBounds, clusteringAssignment, 0));

    // ----------------------------------------------------------------------
    // Make Trees Output
    // ----------------------------------------------------------------------
    printMsg("--- Make Trees Output", ttk::debug::Priority::VERBOSE);
    std::vector<FTMTree_MT *> treesOri(trees);
    if(ShiftMode == 1) { // Star Barycenter
      trees.clear();
      clusteringAssignment.clear();
      for(unsigned int j = 0; j < barycenters.size(); ++j) {
        trees.emplace_back(barycenters[j]);
        clusteringAssignment.emplace_back(j);
      }
      numInputs = trees.size();
    }
    // - Declare VTK arrays
    vtkNew<vtkUnstructuredGrid> vtkArcs{};
    vtkNew<vtkPoints> points{};

    // Node fields
    vtkNew<vtkIntArray> criticalType{};
    criticalType->SetName(ttk::PersistenceCriticalTypeName);
    vtkNew<vtkFloatArray> persistenceNode{};
    persistenceNode->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDNode{};
    clusterIDNode->SetName("ClusterID");
    vtkNew<vtkIntArray> isDummyNode{};
    isDummyNode->SetName("isDummyNode");
    vtkNew<vtkIntArray> branchNodeID{};
    branchNodeID->SetName("BranchNodeID");
    vtkNew<vtkIntArray> pathNodeID{};
    pathNodeID->SetName("PathNodeID");
    vtkNew<vtkFloatArray> scalar{};
    scalar->SetName("Scalar");
    vtkNew<vtkIntArray> isImportantPairsNode{};
    isImportantPairsNode->SetName("isImportantPair");
    vtkNew<vtkIntArray> nodeID{};
    nodeID->SetName("NodeId"); // Simplex Id
    vtkNew<vtkIntArray> trueNodeID{};
    trueNodeID->SetName("TrueNodeId");
    vtkNew<vtkIntArray> vertexID{};
    vertexID->SetName(
      (isPersistenceDiagram ? ttk::VertexScalarFieldName : "VertexId"));

    vtkNew<vtkIntArray> treeIDNode{};
    treeIDNode->SetName("TreeID");
    vtkNew<vtkIntArray> branchBaryNodeID{};
    branchBaryNodeID->SetName("BranchBaryNodeID");
    vtkNew<vtkIntArray> pathBaryNodeID{};
    pathBaryNodeID->SetName("PathBaryNodeID");
    vtkNew<vtkIntArray> isInterpolatedTreeNode{};
    isInterpolatedTreeNode->SetName("isInterpolatedTree");

    vtkNew<vtkFloatArray> percentMatch{};
    percentMatch->SetName("PercentMatchNode");
    vtkNew<vtkFloatArray> persistenceBaryNode{};
    persistenceBaryNode->SetName("PersistenceBarycenter");
    vtkNew<vtkIntArray> persistenceBaryOrderNode{};
    persistenceBaryOrderNode->SetName("PersistenceBarycenterOrder");

    vtkNew<vtkFloatArray> treeNodeId{};
    treeNodeId->SetName("TreeNodeId");

    vtkNew<vtkFloatArray> treeNodeIdOrigin{};
    treeNodeIdOrigin->SetName("TreeNodeIdOrigin");
    vtkNew<vtkDoubleArray> coordinates{};
    coordinates->SetName(ttk::PersistenceCoordinatesName);
    coordinates->SetNumberOfComponents(3);

    vtkNew<vtkIntArray> isMultiPersPairNode{};
    isMultiPersPairNode->SetName("isMultiPersPairNode");

    std::vector<std::vector<double>> customArraysValues(customArrays.size());
    std::vector<std::vector<int>> customIntArraysValues(customIntArrays.size());
    std::vector<std::vector<std::string>> customStringArraysValues(
      customStringArrays.size());

    // Arc fields
    vtkNew<vtkFloatArray> persistenceArc{};
    persistenceArc->SetName("Persistence");
    vtkNew<vtkIntArray> clusterIDArc{};
    clusterIDArc->SetName("ClusterID");
    vtkNew<vtkIntArray> isImportantPairsArc{};
    isImportantPairsArc->SetName("isImportantPair");
    vtkNew<vtkIntArray> isDummyArc{};
    isDummyArc->SetName("isDummyArc");
    vtkNew<vtkIntArray> branchID{};
    branchID->SetName("BranchID");
    vtkNew<vtkIntArray> pathID{};
    pathID->SetName("PathID");
    vtkNew<vtkIntArray> upNodeId{};
    upNodeId->SetName("upNodeId");
    vtkNew<vtkIntArray> downNodeId{};
    downNodeId->SetName("downNodeId");

    vtkNew<vtkIntArray> treeIDArc{};
    treeIDArc->SetName((isPersistenceDiagram ? "DiagramID" : "TreeID"));
    vtkNew<vtkIntArray> branchBaryID{};
    branchBaryID->SetName("BranchBaryNodeID");
    vtkNew<vtkIntArray> pathBaryID{};
    pathBaryID->SetName("PathBaryNodeID");
    vtkNew<vtkIntArray> isInterpolatedTreeArc{};
    isInterpolatedTreeArc->SetName("isInterpolatedTree");

    vtkNew<vtkFloatArray> percentMatchArc{};
    percentMatchArc->SetName("PercentMatchArc");
    vtkNew<vtkFloatArray> persistenceBaryArc{};
    persistenceBaryArc->SetName("PersistenceBarycenter");
    vtkNew<vtkIntArray> persistenceBaryOrderArc{};
    persistenceBaryOrderArc->SetName("PersistenceBarycenterOrder");

    vtkNew<vtkIntArray> pairIdentifier{};
    pairIdentifier->SetName(ttk::PersistencePairIdentifierName);
    vtkNew<vtkIntArray> pairType{};
    pairType->SetName(ttk::PersistencePairTypeName);
    vtkNew<vtkIntArray> pairIsFinite{};
    pairIsFinite->SetName(ttk::PersistenceIsFinite);
    vtkNew<vtkDoubleArray> pairPersistence{};
    pairPersistence->SetName(ttk::PersistenceName);
    vtkNew<vtkDoubleArray> pairBirth{};
    pairBirth->SetName(ttk::PersistenceBirthName);

    vtkNew<vtkIntArray> isMultiPersPairArc{};
    isMultiPersPairArc->SetName("isMultiPersPairArc");

    std::vector<std::vector<double>> customCellArraysValues(
      customArrays.size());
    std::vector<std::vector<int>> customCellIntArraysValues(
      customIntArrays.size());
    std::vector<std::vector<std::string>> customCellStringArraysValues(
      customStringArrays.size());

    // Segmentation
    vtkNew<vtkAppendFilter> appendFilter{};

    // Internal data
    int cellCount = 0;
    int pointCount = 0;
    bool foundOneInterpolatedTree = false;
    nodeCorr.clear();
    nodeCorr.resize(numInputs);
    clusterShift.clear();
    clusterShift.resize(NumberOfBarycenters, 0);
    allBaryPercentMatch.clear();
    allBaryPercentMatch.resize(NumberOfBarycenters);

    // --------------------------------------------------------
    // Iterate through all clusters
    // --------------------------------------------------------
    printMsg("Iterate through all clusters", ttk::debug::Priority::VERBOSE);
    double importantPairsOriginal = importantPairs_;
    for(int c = 0; c < NumberOfBarycenters; ++c) {

      // Get persistence order
      std::vector<int> baryPersistenceOrder;
      if(clusteringOutput and ShiftMode != 1) {
        baryPersistenceOrder.resize(barycenters[c]->getNumberOfNodes(), -1);
        std::vector<std::tuple<ttk::ftm::idNode, ttk::ftm::idNode, dataType>>
          pairsBary;
        barycenters[c]->getPersistencePairsFromTree<dataType>(pairsBary, false);
        for(unsigned int j = 0; j < pairsBary.size(); ++j) {
          int index = pairsBary.size() - 1 - j;
          baryPersistenceOrder[std::get<0>(pairsBary[j])] = index;
          baryPersistenceOrder[std::get<1>(pairsBary[j])] = index;
        }
      }

      // Get radius
      printMsg("// Get radius", ttk::debug::Priority::VERBOSE);
      double delta_max = 1.0;
      int noSample = 0 + noSampleOffset;
      for(int i = 0; i < numInputsOri; ++i) {
        delta_max = std::max(
          (std::get<3>(allBounds[i]) - std::get<2>(allBounds[i])), delta_max);
        delta_max = std::max(
          (std::get<1>(allBounds[i]) - std::get<0>(allBounds[i])), delta_max);
        if(clusteringAssignment[i] != c)
          continue;
        noSample += 1;
      }
      double radius = delta_max * 2 * DimensionSpacing;
      int iSample = 0 + iSampleOffset - 1;

      if(c < NumberOfBarycenters - 1)
        clusterShift[c + 1] = radius * 4 + clusterShift[c];

      // Line/Double line attributes
      prevXMax = 0 + prevXMaxOffset;
      std::vector<double> allPrevXMax;
      double prevYMax = std::numeric_limits<double>::lowest();

      // ------------------------------------------
      // Iterate through all trees of this cluster
      // ------------------------------------------
      printMsg("Iterate through all trees of this cluster",
               ttk::debug::Priority::VERBOSE);
      for(int i = 0; i < numInputs; ++i) {
        if(clusteringAssignment[i] != c)
          continue;

        iSample += 1;

        if((printTreeId == -1 and printClusterId != -1 and c != printClusterId)
           or (printTreeId != -1 and printClusterId == -1 and i != printTreeId)
           or (printTreeId != -1 and printClusterId != -1
               and (c != printClusterId or i != printTreeId)))
          continue;

        // Manage important pairs threshold
        auto fixImportantPairsThreshold
          = [&importantPairsOriginal, this](FTMTree_MT *tree) {
              double importantPairs = importantPairsOriginal;
              if(MaximumImportantPairs > 0 or MinimumImportantPairs > 0) {
                std::vector<std::tuple<idNode, idNode, dataType>> pairs;
                tree->getPersistencePairsFromTree(pairs, false);
                if(MaximumImportantPairs > 0) {
                  int firstIndex = pairs.size() - MaximumImportantPairs;
                  firstIndex
                    = std::max(std::min(firstIndex, int(pairs.size()) - 1), 0);
                  double tempThreshold = 0.999 * std::get<2>(pairs[firstIndex])
                                         / std::get<2>(pairs[pairs.size() - 1]);
                  tempThreshold *= 100;
                  importantPairs = std::max(importantPairs, tempThreshold);
                }
                if(MinimumImportantPairs > 0) {
                  int firstIndex = pairs.size() - MinimumImportantPairs;
                  firstIndex
                    = std::max(std::min(firstIndex, int(pairs.size()) - 1), 0);
                  double tempThreshold = 0.999 * std::get<2>(pairs[firstIndex])
                                         / std::get<2>(pairs[pairs.size() - 1]);
                  tempThreshold *= 100;
                  importantPairs = std::min(importantPairs, tempThreshold);
                }
              }
              return importantPairs;
            };
        importantPairs_ = fixImportantPairsThreshold(trees[i]);

        // Get is interpolated tree (temporal subsampling)
        bool isInterpolatedTree = false;
        if(interpolatedTrees.size() != 0)
          isInterpolatedTree = interpolatedTrees[i];
        foundOneInterpolatedTree |= isInterpolatedTree;

        // Get branching
        printMsg("// Get branching", ttk::debug::Priority::VERBOSE);
        std::vector<idNode> treeBranching, treePathing;
        std::vector<int> treeBranchingID, treePathingID;
        std::vector<bool> isRootPath;
        std::vector<ttk::ftm::idNode> pathOrigin;
        if(not isPersistenceDiagram) {
          trees[i]->getTreeBranching(treeBranching, treeBranchingID);
          if(pathPlanarLayout_ and !pathMatchings.empty()
             and !pathMatchings[c].empty()) {
            isRootPath.resize(trees[i]->getNumberOfNodes());
            std::fill(isRootPath.begin(), isRootPath.end(), false);
            pathOrigin.resize(trees[i]->getNumberOfNodes());
            if(ShiftMode == 1)
              getTreePathing(
                trees[i], pathMatchings[c], true, treePathing, treePathingID);
            else
              getTreePathing(trees[i], pathMatchings[c][i], false, treePathing,
                             treePathingID);
            bool isFirstTree = (ShiftMode == 1);
            std::vector<int> nodeLevel;
            trees[i]->getAllNodeLevel(nodeLevel);
            for(unsigned int j = 0; j < pathMatchings[c].size(); ++j) {
              if((int)j != i and ShiftMode != 1)
                continue;
              for(auto &match : pathMatchings[c][j]) {
                auto first
                  = (isFirstTree ? match.first.first : match.second.first);
                auto second
                  = (isFirstTree ? match.first.second : match.second.second);
                auto lowest
                  = (nodeLevel[first] < nodeLevel[second] ? second : first);
                auto highest
                  = (nodeLevel[first] < nodeLevel[second] ? first : second);
                isRootPath[highest] = true;
                pathOrigin[lowest] = highest;
                pathOrigin[highest] = lowest;
              }
            }
          }
        }

        // Get shift
        printMsg("// Get shift", ttk::debug::Priority::VERBOSE);
        double angle = 360.0 / noSample * iSample;
        double pi = M_PI;
        double diff_x = 0, diff_y = 0;
        double alphaShift
          = BarycenterPositionAlpha ? (-radius + 2 * radius * Alpha) * -1 : 0;
        switch(ShiftMode) {
          case -1:
            diff_x = 0.0;
            diff_y = 0.0;
            break;
          case 0: // Star
            diff_x
              = -1 * radius * std::cos(-1 * angle * pi / 180) + clusterShift[c];
            diff_y = -1 * radius * std::sin(-1 * angle * pi / 180);
            break;
          case 1: // Star Barycenter
            diff_x = clusterShift[c] + alphaShift;
            diff_y = 0;
            break;
          case 2: // Line
            diff_x = prevXMax + radius;
            break;
          case 3: // Double Line
            diff_x = prevXMax + radius;
            if(i >= numInputs / 2) {
              diff_y = -(prevYMax + radius / 2);
              diff_x = allPrevXMax[i - int(numInputs / 2)] + radius;
            } else
              allPrevXMax.emplace_back(prevXMax);
            break;
          default:
            break;
        }

        // isImportantPairVector
        auto getIsImportantPairVector
          = [this](FTMTree_MT *tree, std::vector<bool> &isImportantPairVector,
                   double importantPairs) {
              isImportantPairVector.resize(tree->getNumberOfNodes());
              for(unsigned int n = 0; n < isImportantPairVector.size(); ++n)
                isImportantPairVector[n] = tree->isImportantPair<dataType>(
                  n, importantPairs, excludeImportantPairsLowerValues_,
                  excludeImportantPairsHigherValues_);
            };
        std::vector<bool> isImportantPairVector;
        getIsImportantPairVector(
          trees[i], isImportantPairVector, importantPairs_);

        // Layout correspondence function
        auto getLayoutCorr
          = [](FTMTree_MT *tree, std::vector<SimplexId> &layoutCorr) {
              layoutCorr.resize(tree->getNumberOfNodes());
              int cptNode = 0;
              std::queue<idNode> queueLayoutCorr;
              queueLayoutCorr.emplace(tree->getRoot());
              while(!queueLayoutCorr.empty()) {
                idNode node = queueLayoutCorr.front();
                queueLayoutCorr.pop();

                // Push children to the queue
                std::vector<idNode> children;
                tree->getChildren(node, children);
                for(auto child : children)
                  queueLayoutCorr.emplace(child);

                layoutCorr[node] = cptNode;
                cptNode += 2;
              }
            };

        // Planar layout
        printMsg("// Planar Layout", ttk::debug::Priority::VERBOSE);
        std::vector<SimplexId> layoutCorr;
        getLayoutCorr(trees[i], layoutCorr);
        std::vector<float> layout;
        // TODO remove baryMatchingVector if it is not used
        std::vector<ttk::ftm::idNode> baryMatchingVector;
        // TODO put this declaration just before it is used if the vector is
        // not needed elsewhere
        std::vector<bool> isImportantPairBaryVector;
        if(PlanarLayout) {
          double refPersistence;
          if(clusteringOutput)
            refPersistence = barycenters[0]->getNodePersistence<dataType>(
              barycenters[0]->getRoot());
          else
            refPersistence
              = trees[0]->getNodePersistence<dataType>(trees[0]->getRoot());
          if(not isPersistenceDiagram) {
            treePlanarLayout<dataType>(
              trees[i], allBaryBounds[c], refPersistence, layout);
            // Planar Layout alignment given barycenter
            if(alignTrees and ShiftMode != 1) {
              // Create barycenter layout
              std::vector<float> layoutBary;
              treePlanarLayout<dataType>(
                barycenters[0], allBaryBounds[c], refPersistence, layoutBary);
              std::vector<SimplexId> layoutBaryCorr;
              getLayoutCorr(barycenters[0], layoutBaryCorr);
              // Get barycenter important pairs bool vector
              double isImportantPairBary
                = fixImportantPairsThreshold(barycenters[0]);
              getIsImportantPairVector(
                barycenters[0], isImportantPairBaryVector, isImportantPairBary);
              baryMatchingVector.resize(barycenters[0]->getNumberOfNodes(), -1);
              for(auto match : outputMatchingBarycenter[0][i])
                baryMatchingVector[std::get<0>(match)] = std::get<1>(match);
              std::vector<float> shifts(trees[i]->getNumberOfNodes(), 0);
              for(auto match : outputMatchingBarycenter[0][i]) {
                // Update tree important pair according barycenter
                isImportantPairVector[std::get<1>(match)]
                  = isImportantPairBaryVector[std::get<0>(match)];
                // Update tree layout according barycenter layout
                layout[layoutCorr[std::get<1>(match)]]
                  = layoutBary[layoutBaryCorr[std::get<0>(match)]];
                if(not isImportantPairVector[std::get<1>(match)]) {
                  // Search for birth swap with an important pair
                  ttk::ftm::idNode node = std::get<0>(match);
                  while(not isImportantPairBaryVector[node])
                    node = barycenters[0]->getParentSafe(node);
                  bool baryNodeSup
                    = barycenters[0]->getValue<dataType>(node)
                      > barycenters[0]->getValue<dataType>(std::get<0>(match));
                  bool treeNodeSup
                    = trees[i]->getValue<dataType>(baryMatchingVector[node])
                      > trees[i]->getValue<dataType>(std::get<1>(match));
                  if((not baryNodeSup and treeNodeSup)
                     or (baryNodeSup and not treeNodeSup)) {
                    float shift
                      = std::abs(layout[layoutCorr[std::get<1>(match)]]
                                 - layoutBary[layoutBaryCorr[node]]);
                    layout[layoutCorr[std::get<1>(match)]]
                      = layoutBary[layoutBaryCorr[node]];
                    std::queue<ttk::ftm::idNode> queueShift;
                    queueShift.emplace(
                      trees[i]->getNode(std::get<1>(match))->getOrigin());
                    while(!queueShift.empty()) {
                      ttk::ftm::idNode nodeToShift = queueShift.front();
                      queueShift.pop();
                      shifts[nodeToShift] += shift;
                      ttk::ftm::idNode nodeToShiftParent
                        = trees[i]->getParentSafe(nodeToShift);
                      if(nodeToShiftParent != std::get<1>(match))
                        queueShift.emplace(nodeToShiftParent);
                      std::vector<ttk::ftm::idNode> children;
                      trees[i]->getChildren(nodeToShift, children);
                      for(auto &child : children)
                        queueShift.emplace(child);
                    }
                  }
                }
              }
              for(unsigned int s = 0; s < shifts.size(); ++s)
                layout[layoutCorr[s]] += shifts[s];
            }
          } else {
            persistenceDiagramPlanarLayout<dataType>(trees[i], layout);
          }
        }

        // Get dimension shift
        printMsg("// Get dimension shift", ttk::debug::Priority::VERBOSE);
        double diff_z = PlanarLayout ? 0 : -std::get<4>(allBounds[i]);
        if(DimensionToShift != 0) { // is not X
          float minX = 0;
          if(PlanarLayout) {
            minX = layout[0];
            for(unsigned int l = 0; l < layout.size(); ++l) {
              if(l % 2 == 0)
                minX = std::min(minX, layout[l]);
            }
          }
          double new_diff_x = PlanarLayout ? -minX : -std::get<0>(allBounds[i]);
          bool diffYAllowed
            = (not clusteringOutput
               or (trees.size() == 2 and barycenters.size() == 1));
          if(DimensionToShift == 2) {
            // is Z
            diff_z = -diff_x;
            diff_x = new_diff_x;
          } else if(diffYAllowed and DimensionToShift == 1) {
            // is Y
            diff_y = diff_x;
            diff_x = new_diff_x;
          } else if(DimensionToShift == 3) {
            // Custom
            if(diffYAllowed)
              diff_y = YShift * diff_x + (1 - YShift) * diff_y;
            diff_z = ZShift * -diff_x + (1 - ZShift) * diff_z;
            diff_x = XShift * diff_x + (1 - XShift) * new_diff_x;
          }
        }

        // Internal arrays
        printMsg("// Internal arrays", ttk::debug::Priority::VERBOSE);
        nodeCorr[i].resize(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> treeSimplexId(trees[i]->getNumberOfNodes());
        std::vector<SimplexId> treeDummySimplexId(trees[i]->getNumberOfNodes());
        std::vector<idNode> treeMatching(trees[i]->getNumberOfNodes(), -1);
        if(clusteringOutput and ShiftMode != 1)
          for(auto match : outputMatchingBarycenter[c][i])
            treeMatching[std::get<1>(match)] = std::get<0>(match);
        // _ m[i][j] contains the node in treesOri[j] matched to the node i in
        // the barycenter
        std::vector<std::vector<idNode>> baryMatching(
          trees[i]->getNumberOfNodes(), std::vector<idNode>(numInputsOri, -1));
        if(ShiftMode == 1) {
          for(size_t j = 0; j < outputMatchingBarycenter[c].size(); ++j)
            for(auto match : outputMatchingBarycenter[c][j])
              baryMatching[std::get<0>(match)][j] = std::get<1>(match);
          allBaryPercentMatch[c].resize(trees[i]->getNumberOfNodes(), 100.0);
        }
        double minBirth = std::numeric_limits<double>::max(),
               maxBirth = std::numeric_limits<double>::lowest();
        SimplexId minBirthNode = 0, maxBirthNode = 0;

        // ----------------------------
        // Tree traversal
        // ----------------------------
        printMsg("// Tree traversal", ttk::debug::Priority::VERBOSE);
        std::queue<idNode> queue;
        queue.emplace(trees[i]->getRoot());
        while(!queue.empty()) {
          idNode node = queue.front();
          queue.pop();
          idNode nodeOrigin = trees[i]->getNode(node)->getOrigin();
          idNode nodeParent = trees[i]->getParentSafe(node);

          // Push children to the queue
          printMsg(
            "// Push children to the queue", ttk::debug::Priority::VERBOSE);
          std::vector<idNode> children;
          trees[i]->getChildren(node, children);
          for(auto child : children)
            queue.emplace(child);

          // --------------
          // Insert point
          // --------------
          auto getPoint
            = [&](vtkUnstructuredGrid *vtu, int pointID, double(&point)[3]) {
                if(not isPersistenceDiagram) {
                  double *pointTemp = vtu->GetPoints()->GetPoint(pointID);
                  for(int k = 0; k < 3; ++k)
                    point[k] += pointTemp[k];
                } else {
                  for(int k = 0; k < 3; ++k)
                    point[k] = vtu->GetPointData()
                                 ->GetArray(ttk::PersistenceCoordinatesName)
                                 ->GetComponent(pointID, k);
                }
              };

          printMsg("// Get and insert point", ttk::debug::Priority::VERBOSE);
          int nodeMesh = -1;
          int nodeMeshTreeIndex = -1;
          double noMatched = 0.0;
          double point[3] = {0, 0, 0};
          if(ShiftMode == 1) { // Star barycenter
            for(int j = 0; j < numInputsOri; ++j) {
              if((int)baryMatching[node][j] != -1) {
                nodeMesh = treesNodeCorrMesh[j][baryMatching[node][j]];
                getPoint(treesNodes[j], nodeMesh, point);
                noMatched += 1;
                nodeMeshTreeIndex = j;
              }
            }
            for(int k = 0; k < 3; ++k)
              point[k] /= noMatched;
          } else if(not isInterpolatedTree and treesNodes.size() != 0
                    and treesNodes[i] != nullptr) {
            nodeMesh = treesNodeCorrMesh[i][node];
            getPoint(treesNodes[i], nodeMesh, point);
          }
          if(PlanarLayout) {
            point[0] = layout[layoutCorr[node]];
            point[1] = layout[layoutCorr[node] + 1];
            point[2] = 0;
          }
          point[0] += diff_x;
          point[1] += diff_y;
          point[2] += diff_z;

          // Bary percentage matching
          if(ShiftMode == 1) { // Star Barycenter
            float percentMatchT = noMatched * 100 / numInputs;
            allBaryPercentMatch[c][node] = percentMatchT;
          }

          // Get x Max and y Min for next iteration if needed (double line mode)
          prevXMax = std::max(prevXMax, point[0]);
          if(ShiftMode == 3) { // Double line
            if(i < numInputs / 2)
              prevYMax = std::max(prevYMax, point[1]);
            if(i == int(numInputs / 2) - 1)
              prevXMax = 0;
          }

          bool dummyNode
            = PlanarLayout and not branchDecompositionPlanarLayout_
              and ((!trees[i]->isRoot(node) and !trees[i]->isLeaf(node))
                   or isPersistenceDiagram);
          if(dummyNode) {
            double pointToAdd[3];
            if(not isPersistenceDiagram)
              // will be modified when processing son
              std::copy(
                std::begin(point), std::end(point), std::begin(pointToAdd));
            else {
              double pdPoint[3] = {
                point[0],
                point[1]
                  - (layout[layoutCorr[node] + 1] - layout[layoutCorr[node]]),
                0};
              std::copy(
                std::begin(pdPoint), std::end(pdPoint), std::begin(pointToAdd));
            }
            treeDummySimplexId[node] = points->InsertNextPoint(pointToAdd);
            if(isPersistenceDiagram) {
              if(layout[layoutCorr[node]] < minBirth) {
                minBirth = layout[layoutCorr[node]];
                minBirthNode = treeDummySimplexId[node];
              }
              if(layout[layoutCorr[node]] > maxBirth) {
                maxBirth = layout[layoutCorr[node]];
                maxBirthNode = treeDummySimplexId[node];
              }
            }
          }
          SimplexId nextPointId = points->InsertNextPoint(point);
          treeSimplexId[node] = nextPointId;
          nodeCorr[i][node] = nextPointId;
          if(dummyNode and not pathPlanarLayout_)
            nodeCorr[i][node] = treeDummySimplexId[node];
          if(isPersistenceDiagram)
            nodeCorr[i][node] = nextPointId;

          idNode nodeBranching
            = ((PlanarLayout and branchDecompositionPlanarLayout_)
                   or isPersistenceDiagram
                 ? node
                 : treeBranching[node]);

          // Path Layout Dummy Node
          bool pathDummyNode = false;

          // --------------
          // Insert cell connecting parent
          // --------------
          printMsg(
            "// Add cell connecting parent", ttk::debug::Priority::VERBOSE);
          if(!trees[i]->isRoot(node) or isPersistenceDiagram) {
            vtkIdType pointIds[2];
            pointIds[0] = treeSimplexId[node];

            // TODO too many dummy cells are created
            bool dummyCell
              = PlanarLayout and not branchDecompositionPlanarLayout_
                and (node < treeBranching.size()
                     and treeBranching[node] == nodeParent)
                and !trees[i]->isRoot(nodeParent) and not isPersistenceDiagram;
            if(isPersistenceDiagram) {
              pointIds[1] = treeDummySimplexId[node];
            } else if(PlanarLayout and branchDecompositionPlanarLayout_) {
              pointIds[1] = treeSimplexId[treeBranching[node]];
            } else if(dummyCell) {
              double dummyPoint[3]
                = {point[0], layout[layoutCorr[nodeParent] + 1] + diff_y,
                   0. + diff_z};
              SimplexId dummyPointId = treeDummySimplexId[nodeParent];
              points->SetPoint(dummyPointId, dummyPoint);
              vtkIdType dummyPointIds[2];
              dummyPointIds[0] = dummyPointId;
              dummyPointIds[1] = treeSimplexId[nodeParent];
              vtkArcs->InsertNextCell(VTK_LINE, 2, dummyPointIds);
              pointIds[1] = dummyPointId;
            } else
              pointIds[1] = treeSimplexId[nodeParent];

            // Path Layout Dummy Cell
            bool isNodeParentImportant = isImportantPairVector[nodeParent];
            bool pathDummyCell = not dummyCell and pathPlanarLayout_
                                 and isNodeParentImportant
                                 and !trees[i]->isRoot(nodeParent);
            if(not pathDummyCell and alignTrees and ShiftMode != 1) {
              pathDummyCell
                = not dummyCell and pathPlanarLayout_
                  and layout[layoutCorr[node]] != layout[layoutCorr[nodeParent]]
                  and !trees[i]->isRoot(nodeParent);
            }
            if(pathDummyCell) {
              pathDummyNode = true;
              double pathDummyPoint[3]
                = {layout[layoutCorr[node]] + diff_x,
                   layout[layoutCorr[nodeParent] + 1] + diff_y, 0. + diff_z};
              SimplexId pathDummyPointId
                = points->InsertNextPoint(pathDummyPoint);
              vtkIdType pathDummyCellPointIds[2];
              pathDummyCellPointIds[0] = pathDummyPointId;
              pathDummyCellPointIds[1] = treeSimplexId[nodeParent];
              vtkArcs->InsertNextCell(VTK_LINE, 2, pathDummyCellPointIds);
              pointIds[1] = pathDummyPointId;
            }
            vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);

            // --------------
            // Arc field
            // --------------
            int toAdd = 1 + (dummyCell ? 1 : 0) + (pathDummyCell ? 1 : 0);
            for(int toAddT = 0; toAddT < toAdd; ++toAddT) {
              // Add arc matching percentage
              if(ShiftMode == 1) { // Star Barycenter
                auto nodeToGet
                  = (!isPersistenceDiagram ? allBaryBranching[c][node] : node);
                percentMatchArc->InsertNextTuple1(
                  allBaryPercentMatch[c][nodeToGet]);
              }

              // Add branch/path bary ID
              printMsg(
                "// Push arc bary branch id", ttk::debug::Priority::VERBOSE);
              if(clusteringOutput and ShiftMode != 1) {
                // Branch
                int tBranchID = -1;
                auto nodeToGet = node;
                if(treeMatching[nodeToGet] < allBaryBranchingID[c].size())
                  tBranchID = allBaryBranchingID[c][treeMatching[nodeToGet]];
                branchBaryID->InsertNextTuple1(tBranchID);
                // Path
                if(pathPlanarLayout_ and !pathMatchings.empty()
                   and !pathMatchings[c].empty()) {
                  int tPathID = -1;
                  nodeToGet = node;
                  if(treeMatching[nodeToGet] < allBaryPathingID[c].size())
                    tPathID = allBaryPathingID[c][treeMatching[nodeToGet]];
                  pathBaryID->InsertNextTuple1(tPathID);
                }
              }

              // Add branch/path ID
              if(not isPersistenceDiagram) {
                // Branch
                int tBranchID = treeBranchingID[node];
                branchID->InsertNextTuple1(tBranchID);
                // Path
                if(pathPlanarLayout_ and !pathMatchings.empty()
                   and !pathMatchings[c].empty()) {
                  int tPathID = treePathingID[node];
                  pathID->InsertNextTuple1(tPathID);
                }
              }

              // Add up and down nodeId
              if(not isPersistenceDiagram) {
                upNodeId->InsertNextTuple1(treeSimplexId[nodeParent]);
                downNodeId->InsertNextTuple1(treeSimplexId[node]);
              }

              // Add arc persistence
              printMsg(
                "// Push arc persistence", ttk::debug::Priority::VERBOSE);
              idNode nodeToGetPers = nodeBranching;
              double persToAdd
                = trees[i]->getNodePersistence<dataType>(nodeToGetPers);
              persistenceArc->InsertNextTuple1(persToAdd);

              // Add arc persistence barycenter and order
              if(clusteringOutput and ShiftMode != 1) {
                idNode nodeToGet = nodeBranching;
                if(treeMatching[nodeToGet] < allBaryBranchingID[c].size()) {
                  persistenceBaryArc->InsertTuple1(
                    cellCount, barycenters[c]->getNodePersistence<dataType>(
                                 treeMatching[nodeToGet]));
                  persistenceBaryOrderArc->InsertTuple1(
                    cellCount, baryPersistenceOrder[treeMatching[nodeToGet]]);
                }
              }

              // Add arc cluster ID
              clusterIDArc->InsertNextTuple1(clusteringAssignment[i]);

              // Add arc tree ID
              treeIDArc->InsertNextTuple1(i);

              // Add isImportantPair
              bool isImportant = isImportantPairVector[nodeBranching];
              isImportantPairsArc->InsertNextTuple1(isImportant);

              // Add isDummyArc
              bool isDummy = toAdd >= 2 and toAddT == (toAdd - 2);
              isDummyArc->InsertNextTuple1(isDummy);

              // Add isInterpolatedTree
              isInterpolatedTreeArc->InsertNextTuple1(isInterpolatedTree);

              // Add pairIdentifier
              pairIdentifier->InsertNextTuple1(treeSimplexId[node]);

              // Add birth and death
              auto birthDeath = trees[i]->getBirthDeath<dataType>(node);
              pairBirth->InsertNextTuple1(std::get<0>(birthDeath));
              pairPersistence->InsertNextTuple1(std::get<1>(birthDeath)
                                                - std::get<0>(birthDeath));

              // Add isMinMaxPair
              bool isMinMaxPair
                = (trees[i]->isRoot(node) and not trees[i]->isLeaf(node))
                  or (trees[i]->isRoot(nodeOrigin)
                      and not trees[i]->isLeaf(nodeOrigin));
              pairIsFinite->InsertNextTuple1(!isMinMaxPair);

              // Add pairType TODO
              pairType->InsertNextTuple1(0);

              // Add isMultiPersPairArc
              bool isMultiPersPair
                = (trees[i]->isMultiPersPair(nodeBranching)
                   or trees[i]->isMultiPersPair(
                     trees[i]->getNode(nodeBranching)->getOrigin()));
              isMultiPersPairArc->InsertNextTuple1(isMultiPersPair);

              // Add custom point arrays to cells
              for(unsigned int ca = 0; ca < customArrays.size(); ++ca)
                customCellArraysValues[ca].push_back(
                  std::get<1>(customArrays[ca])[nodeBranching]);
              for(unsigned int ca = 0; ca < customIntArrays.size(); ++ca)
                customCellIntArraysValues[ca].push_back(
                  std::get<1>(customIntArrays[ca])[nodeBranching]);
              for(unsigned int ca = 0; ca < customStringArrays.size(); ++ca)
                customCellStringArraysValues[ca].push_back(
                  std::get<1>(customStringArrays[ca])[nodeBranching]);

              cellCount++;
            }
          }

          // --------------
          // Node field
          // --------------
          int toAdd = 1 + (dummyNode ? 1 : 0) + (pathDummyNode ? 1 : 0);
          for(int toAddT = 0; toAddT < toAdd; ++toAddT) {
            bool isPathDummyNode = pathDummyNode and toAdd >= 2
                                   and toAddT == (toAdd - 1)
                                   and !trees[i]->isRoot(node);
            auto nodeToGet = (isPathDummyNode ? nodeParent : node);
            // Add node id
            nodeID->InsertNextTuple1(treeSimplexId[nodeToGet]);

            // Add trueNodeId
            trueNodeID->InsertNextTuple1(nodeToGet);

            // Add VertexId
            int nodeVertexId = -1;
            if(i < int(treesNodes.size()) and treesNodes[i]) {
              auto vertexIdArray = treesNodes[i]->GetPointData()->GetArray(
                (isPersistenceDiagram ? ttk::VertexScalarFieldName
                                      : "VertexId"));
              if(vertexIdArray and nodeMesh != -1)
                nodeVertexId = vertexIdArray->GetTuple1(nodeMesh);
            }
            vertexID->InsertNextTuple1(nodeVertexId);

            // Add node scalar
            auto scalarValue = trees[i]->getValue<dataType>(nodeToGet);
            scalar->InsertNextTuple1(scalarValue);

            // Add criticalType
            printMsg("// Add criticalType", ttk::debug::Priority::VERBOSE);
            int criticalTypeT = -1;
            if(not isPersistenceDiagram) {
              if(not isInterpolatedTree) {
                if(ShiftMode == 1) {
                  if(nodeMeshTreeIndex != -1) {
                    auto array
                      = treesNodes[nodeMeshTreeIndex]->GetPointData()->GetArray(
                        "CriticalType");
                    criticalTypeT = array->GetTuple1(nodeMesh);
                  }
                } else if(treesNodes.size() != 0 and treesNodes[i] != nullptr) {
                  auto array
                    = treesNodes[i]->GetPointData()->GetArray("CriticalType");
                  criticalTypeT = array->GetTuple1(nodeMesh);
                }
              } else {
                // TODO critical type for interpolated trees
              }
            } else {
              auto locMin = static_cast<int>(ttk::CriticalType::Local_minimum);
              auto saddle1 = static_cast<int>(ttk::CriticalType::Saddle1);
              auto locMax = static_cast<int>(ttk::CriticalType::Local_maximum);
              auto saddle2 = static_cast<int>(ttk::CriticalType::Saddle2);
              auto nodeIsRoot = trees[i]->isRoot(node);
              criticalTypeT
                = (toAddT == 1
                     ? (isPDSadMax or nodeIsRoot ? locMax : saddle1)
                     : (not isPDSadMax or nodeIsRoot ? locMin : saddle2));
            }
            criticalType->InsertNextTuple1(criticalTypeT);

            // Add node matching percentage
            if(ShiftMode == 1) // Star Barycenter
              percentMatch->InsertNextTuple1(allBaryPercentMatch[c][node]);

            // Add node branch/path bary id
            printMsg(
              "// Add node bary branch id", ttk::debug::Priority::VERBOSE);
            if(clusteringOutput and ShiftMode != 1) {
              // Branch
              int tBranchID = -1;
              auto branchNode
                = (isPathDummyNode
                     ? trees[i]->getNode(treeBranching[node])->getOrigin()
                     : (trees[i]->isLeaf(node) ? node : nodeOrigin));
              if(treeMatching[branchNode] < allBaryBranchingID[c].size())
                tBranchID = allBaryBranchingID[c][treeMatching[branchNode]];
              branchBaryNodeID->InsertNextTuple1(tBranchID);
              // Path
              if(pathPlanarLayout_ and !pathMatchings.empty()
                 and !pathMatchings[c].empty()) {
                int tPathID = -1;
                auto pathNode
                  = (isPathDummyNode
                       ? node
                       : (not isRootPath[node] ? node : pathOrigin[node]));
                if(treeMatching[pathNode] < allBaryPathingID[c].size())
                  tPathID = allBaryPathingID[c][treeMatching[pathNode]];
                pathBaryNodeID->InsertNextTuple1(tPathID);
              }
            }

            // Add node branch/path id
            if(not isPersistenceDiagram) {
              // Branch
              auto branchNode
                = (isPathDummyNode
                     ? trees[i]->getNode(treeBranching[node])->getOrigin()
                     : (trees[i]->isLeaf(node) ? node : nodeOrigin));
              int tBranchID = treeBranchingID[branchNode];
              branchNodeID->InsertNextTuple1(tBranchID);
              // Path
              if(pathPlanarLayout_ and !pathMatchings.empty()
                 and !pathMatchings[c].empty()) {
                auto pathNode
                  = (isPathDummyNode
                       ? node
                       : (not isRootPath[node] ? node : pathOrigin[node]));
                // auto pathNode = node;
                int tPathID = treePathingID[pathNode];
                pathNodeID->InsertNextTuple1(tPathID);
              }
            }

            // Add node persistence
            printMsg("// Push node persistence", ttk::debug::Priority::VERBOSE);
            persistenceNode->InsertNextTuple1(
              trees[i]->getNodePersistence<dataType>(node));

            // Add node persistence barycenter
            if(clusteringOutput and ShiftMode != 1) {
              if(treeMatching[node] < allBaryBranchingID[c].size()) {
                persistenceBaryNode->InsertTuple1(
                  pointCount, barycenters[c]->getNodePersistence<dataType>(
                                treeMatching[node]));
                persistenceBaryOrderNode->InsertTuple1(
                  pointCount, baryPersistenceOrder[treeMatching[node]]);
              }
            }

            // Add node clusterID
            clusterIDNode->InsertNextTuple1(clusteringAssignment[i]);

            // Add node tree ID
            treeIDNode->InsertNextTuple1(i);

            // Add isDummyNode
            bool isDummy
              = toAdd == 2 and toAddT == 1 and !trees[i]->isRoot(node);
            if(pathPlanarLayout_) {
              isDummy = (toAdd >= 2 and dummyNode and toAddT == 0
                         and !trees[i]->isRoot(node));
              isDummy = isDummy or isPathDummyNode;
            }
            isDummyNode->InsertNextTuple1(isDummy);

            // Add isInterpolatedTree
            isInterpolatedTreeNode->InsertNextTuple1(isInterpolatedTree);

            // Add isImportantPair
            bool isImportant = isImportantPairVector[nodeToGet];
            isImportantPairsNode->InsertNextTuple1(isImportant);

            // Add treeNodeId
            treeNodeId->InsertNextTuple1(nodeToGet);

            // Add treeNodeIdOrigin
            treeNodeIdOrigin->InsertNextTuple1(nodeOrigin);

            // Add coordinates
            if(isPersistenceDiagram and !treesNodes.empty()
               and ShiftMode != 1) {
              double coord[3] = {0.0, 0.0, 0.0};
              getPoint(treesNodes[i], treesNodeCorrMesh[i][node], coord);
              coordinates->InsertNextTuple3(coord[0], coord[1], coord[2]);
            }

            // Add isMultiPersPairArc
            bool isMultiPersPair = (trees[i]->isMultiPersPair(node)
                                    or trees[i]->isMultiPersPair(
                                      trees[i]->getNode(node)->getOrigin()));
            isMultiPersPairNode->InsertNextTuple1(isMultiPersPair);

            // Add custom arrays
            for(unsigned int ca = 0; ca < customArrays.size(); ++ca)
              customArraysValues[ca].emplace_back(
                std::get<1>(customArrays[ca])[nodeToGet]);
            for(unsigned int ca = 0; ca < customIntArrays.size(); ++ca)
              customIntArraysValues[ca].emplace_back(
                std::get<1>(customIntArrays[ca])[nodeToGet]);
            for(unsigned int ca = 0; ca < customStringArrays.size(); ++ca)
              customStringArraysValues[ca].emplace_back(
                std::get<1>(customStringArrays[ca])[nodeToGet]);

            pointCount++;
          }

          printMsg("end loop", ttk::debug::Priority::VERBOSE);
        } // end tree traversal

        // Add diagonal if isPersistenceDiagram
        if(isPersistenceDiagram) {
          vtkIdType pointIds[2];
          pointIds[0] = minBirthNode;
          pointIds[1] = maxBirthNode;
          vtkArcs->InsertNextCell(VTK_LINE, 2, pointIds);
          cellCount++;

          pairIdentifier->InsertNextTuple1(-1);
          pairType->InsertNextTuple1(-1);
          pairPersistence->InsertNextTuple1(-1);
          pairIsFinite->InsertNextTuple1(0);
          pairBirth->InsertNextTuple1(0);

          for(unsigned int ca = 0; ca < customArrays.size(); ++ca)
            customCellArraysValues[ca].push_back(-1);
          for(unsigned int ca = 0; ca < customIntArrays.size(); ++ca)
            customCellIntArraysValues[ca].push_back(-1);
          for(unsigned int ca = 0; ca < customStringArrays.size(); ++ca)
            customCellStringArraysValues[ca].push_back("");

          isMultiPersPairArc->InsertNextTuple1(0);
          clusterIDArc->InsertNextTuple1(clusteringAssignment[i]);
          treeIDArc->InsertNextTuple1(i);
          isImportantPairsArc->InsertNextTuple1(0);
          branchBaryID->InsertNextTuple1(-1);
          percentMatchArc->InsertNextTuple1(100);
        }

        // --------------
        // Manage segmentation
        // --------------
        // Use TransformFilter (see commit
        // 85600763a8907674b8e57d6ad77ca97640725b30) when issue #513 is
        // solved.
        printMsg("// Shift segmentation", ttk::debug::Priority::VERBOSE);
        if(OutputSegmentation and not PlanarLayout and treesSegmentation[i]) {
          vtkNew<vtkUnstructuredGrid> iTreesSegmentationCopy{};
          if(ShiftMode != -1)
            iTreesSegmentationCopy->DeepCopy(treesSegmentation[i]);
          else
            iTreesSegmentationCopy->ShallowCopy(treesSegmentation[i]);
          auto iVkOutputSegmentationTemp
            = vtkUnstructuredGrid::SafeDownCast(iTreesSegmentationCopy);
          if(!iVkOutputSegmentationTemp
             or !iVkOutputSegmentationTemp->GetPoints()) {
            printWrn("Convert segmentation to vtkUnstructuredGrid.");
            vtkNew<vtkAppendFilter> appendFilter2{};
            appendFilter2->AddInputData(treesSegmentation[i]);
            appendFilter2->Update();
            iVkOutputSegmentationTemp->ShallowCopy(appendFilter2->GetOutput());
          }
          if(ShiftMode != -1) {
            for(int p = 0;
                p < iVkOutputSegmentationTemp->GetPoints()->GetNumberOfPoints();
                ++p) {
              double *point
                = iVkOutputSegmentationTemp->GetPoints()->GetPoint(p);
              point[0] += diff_x;
              point[1] += diff_y;
              point[2] += diff_z;
              iVkOutputSegmentationTemp->GetPoints()->SetPoint(p, point);
            }
          }
          appendFilter->AddInputData(iVkOutputSegmentationTemp);
        }
        printMsg("// Shift segmentation DONE", ttk::debug::Priority::VERBOSE);
      }
    }
    for(int i = persistenceBaryNode->GetNumberOfTuples(); i < pointCount; ++i)
      persistenceBaryNode->InsertNextTuple1(0);
    for(int i = persistenceBaryArc->GetNumberOfTuples(); i < cellCount; ++i)
      persistenceBaryArc->InsertNextTuple1(0);
    for(int i = persistenceBaryOrderNode->GetNumberOfTuples(); i < pointCount;
        ++i)
      persistenceBaryOrderNode->InsertNextTuple1(0);
    for(int i = persistenceBaryOrderArc->GetNumberOfTuples(); i < cellCount;
        ++i)
      persistenceBaryOrderArc->InsertNextTuple1(0);

    // --- Add VTK arrays to output
    printMsg("// Add VTK arrays to output", ttk::debug::Priority::VERBOSE);
    // - Manage node output
    // Custom arrays
    addVtkCustomArrays(customArrays, customArraysValues, vtkOutputNode, 0, 0);
    addVtkCustomArrays(
      customIntArrays, customIntArraysValues, vtkOutputNode, 1, 0);
    addVtkCustomArrays(
      customStringArrays, customStringArraysValues, vtkOutputNode, 2, 0);

    // Classical arrays
    vtkOutputNode->SetPoints(points);
    vtkOutputNode->GetPointData()->AddArray(criticalType);
    vtkOutputNode->GetPointData()->AddArray(persistenceNode);
    vtkOutputNode->GetPointData()->AddArray(clusterIDNode);
    vtkOutputNode->GetPointData()->AddArray(treeIDNode);
    vtkOutputNode->GetPointData()->AddArray(trueNodeID);
    vtkOutputNode->GetPointData()->AddArray(vertexID);
    vtkOutputNode->GetPointData()->AddArray(isImportantPairsNode);
    vtkOutputNode->GetPointData()->AddArray(isMultiPersPairNode);
    if(not isPersistenceDiagram) {
      vtkOutputNode->GetPointData()->AddArray(nodeID);
      vtkOutputNode->GetPointData()->AddArray(branchNodeID);
      vtkOutputNode->GetPointData()->AddArray(isDummyNode);
      if(pathPlanarLayout_ and !pathMatchings.empty()
         and !pathMatchings[0].empty())
        vtkOutputNode->GetPointData()->AddArray(pathNodeID);
    }
    if(not branchDecompositionPlanarLayout_ and not isPersistenceDiagram)
      vtkOutputNode->GetPointData()->AddArray(scalar);
    if(clusteringOutput and ShiftMode != 1) {
      vtkOutputNode->GetPointData()->AddArray(branchBaryNodeID);
      if(pathPlanarLayout_ and !pathMatchings.empty()
         and !pathMatchings[0].empty())
        vtkOutputNode->GetPointData()->AddArray(pathBaryNodeID);
      vtkOutputNode->GetPointData()->AddArray(persistenceBaryNode);
      vtkOutputNode->GetPointData()->AddArray(persistenceBaryOrderNode);
    }
    if(foundOneInterpolatedTree)
      vtkOutputNode->GetPointData()->AddArray(isInterpolatedTreeNode);
    if(ShiftMode == 1) // Star Barycenter
      vtkOutputNode->GetPointData()->AddArray(percentMatch);
    if(outputTreeNodeIndex)
      vtkOutputNode->GetPointData()->AddArray(treeNodeId);
    if(isPersistenceDiagram) {
      vtkOutputNode->GetPointData()->AddArray(treeNodeIdOrigin);
      if(!treesNodes.empty() and ShiftMode != 1)
        vtkOutputNode->GetPointData()->AddArray(coordinates);
    }

    // - Manage arc output
    // Custom arrays
    addVtkCustomArrays(customArrays, customCellArraysValues, vtkArcs, 0, 1);
    addVtkCustomArrays(
      customIntArrays, customCellIntArraysValues, vtkArcs, 1, 1);
    addVtkCustomArrays(
      customStringArrays, customCellStringArraysValues, vtkArcs, 2, 1);

    // Classical arrays
    vtkArcs->SetPoints(points);
    vtkArcs->GetCellData()->AddArray(persistenceArc);
    vtkArcs->GetCellData()->AddArray(clusterIDArc);
    vtkArcs->GetCellData()->AddArray(treeIDArc);
    vtkArcs->GetCellData()->AddArray(isImportantPairsArc);
    vtkArcs->GetCellData()->AddArray(isMultiPersPairArc);
    if(not isPersistenceDiagram) {
      vtkArcs->GetCellData()->AddArray(isDummyArc);
      vtkArcs->GetCellData()->AddArray(branchID);
      if(pathPlanarLayout_ and !pathMatchings.empty()
         and !pathMatchings[0].empty())
        vtkArcs->GetCellData()->AddArray(pathID);
      vtkArcs->GetCellData()->AddArray(upNodeId);
      vtkArcs->GetCellData()->AddArray(downNodeId);
    }
    if(clusteringOutput and ShiftMode != 1) {
      vtkArcs->GetCellData()->AddArray(branchBaryID);
      if(pathPlanarLayout_ and !pathMatchings.empty()
         and !pathMatchings[0].empty())
        vtkArcs->GetCellData()->AddArray(pathBaryID);
      vtkArcs->GetCellData()->AddArray(persistenceBaryArc);
      vtkArcs->GetCellData()->AddArray(persistenceBaryOrderArc);
    }
    if(foundOneInterpolatedTree)
      vtkArcs->GetCellData()->AddArray(isInterpolatedTreeArc);
    if(ShiftMode == 1) // Star Barycenter
      vtkArcs->GetCellData()->AddArray(percentMatchArc);
    if(isPersistenceDiagram) {
      vtkArcs->GetCellData()->AddArray(pairIdentifier);
      vtkArcs->GetCellData()->AddArray(pairType);
      vtkArcs->GetCellData()->AddArray(pairIsFinite);
      vtkArcs->GetCellData()->AddArray(pairPersistence);
      vtkArcs->GetCellData()->AddArray(pairBirth);
    }
    if(vtkOutputArc == vtkOutputNode)
      vtkArcs->GetPointData()->ShallowCopy(vtkOutputNode->GetPointData());
    if(not branchDecompositionPlanarLayout_)
      vtkArcs->GetPointData()->AddArray(scalar);
    vtkOutputArc->ShallowCopy(vtkArcs);

    // - Manage segmentation output
    if(OutputSegmentation and not PlanarLayout
       and appendFilter->GetNumberOfInputConnections(0) != 0) {
      appendFilter->SetMergePoints(false);
      appendFilter->Update();
      vtkOutputSegmentation->ShallowCopy(appendFilter->GetOutput());
    }

    //
    if(ShiftMode == 1) // Star Barycenter
      trees = treesOri;
  }

  // ==========================================================================
  // Bounds Utils
  // ==========================================================================
  std::tuple<double, double, double, double, double, double>
    getRealBounds(vtkUnstructuredGrid *treeNodes,
                  FTMTree_MT *tree,
                  std::vector<int> &nodeCorrT) {
    double x_min = std::numeric_limits<double>::max();
    double y_min = std::numeric_limits<double>::max();
    double z_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_max = std::numeric_limits<double>::lowest();
    double z_max = std::numeric_limits<double>::lowest();
    std::queue<idNode> queue;
    queue.emplace(tree->getRoot());
    while(!queue.empty()) {
      idNode node = queue.front();
      queue.pop();
      double *point = treeNodes->GetPoints()->GetPoint(nodeCorrT[node]);
      x_min = std::min(x_min, point[0]);
      x_max = std::max(x_max, point[0]);
      y_min = std::min(y_min, point[1]);
      y_max = std::max(y_max, point[1]);
      z_min = std::min(z_min, point[2]);
      z_max = std::max(z_max, point[2]);
      std::vector<idNode> children;
      tree->getChildren(node, children);
      for(auto child : children)
        queue.emplace(child);
    }
    return std::make_tuple(x_min, x_max, y_min, y_max, z_min, z_max);
  }

  std::tuple<double, double, double, double, double, double>
    getRealBounds(vtkUnstructuredGrid *treeNodes, FTMTree_MT *tree) {
    std::vector<int> nodeCorrT(tree->getNumberOfNodes());
    for(size_t i = 0; i < nodeCorrT.size(); ++i)
      nodeCorrT[i] = i;

    return getRealBounds(treeNodes, tree, nodeCorrT);
  }
};
