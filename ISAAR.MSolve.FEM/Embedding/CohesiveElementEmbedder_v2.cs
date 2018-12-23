using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Collections.Generic;
using System.Linq;
using IEmbeddedElement = ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;


namespace ISAAR.MSolve.FEM.Embedding
{

    /// <summary>
    /// This class should only be used with <see cref="FEM.Elements.cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlk"/>
    /// </summary>
    public class CohesiveElementEmbedder_v2 : IElementDOFEnumerator
    {
        private readonly Model_v2 model;
        private readonly IElement embeddedElement;
        private readonly IEmbeddedDOFInHostTransformationVector transformation;
        private readonly Dictionary<SuperElementDOF, int> superElementMap = new Dictionary<SuperElementDOF, int>();
        private readonly Dictionary<EmbeddedNode, Dictionary<DOFType, int>> dofToHostMapping = new Dictionary<EmbeddedNode, Dictionary<DOFType, int>>();
        private Matrix2D transformationMatrix; //TODO: use sparse CSC matrix for this

        public CohesiveElementEmbedder_v2(Model_v2 model, Element embeddedElement, IEmbeddedDOFInHostTransformationVector transformation)
        {
            this.model = model;
            this.embeddedElement = embeddedElement;
            this.transformation = transformation;
            Initialize();            
        }

       
        private void InitializeMappings()
        {
            var e = embeddedElement.IElementType as IEmbeddedElement;
            superElementMap.Clear();
            int index = 0;
            foreach (var embeddedNode in e.EmbeddedNodes)
            {
                int nodeOrderInEmbeddedElement = embeddedElement.INodes.IndexOf(embeddedNode.Node);
                var currentEmbeddedNodeDOFs = embeddedElement.IElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement];
                //var currentNodeDOFs = currentEmbeddedNodeDOFs.Intersect(embeddedNode.DependentDOFs);
                var independentEmbeddedDOFs = currentEmbeddedNodeDOFs.Except(embeddedNode.DependentDOFs);

                // TODO: Optimization to exclude host DOFs that embedded node does not depend on.
                for (int i = 0; i < embeddedNode.EmbeddedInElement.Nodes.Count; i++)
                {
                    var currentNodeDOFs = embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedNode.EmbeddedInElement)[i];
                    foreach (var dof in currentNodeDOFs)
                    {
                        var superElementDOF = new SuperElementDOF() { DOF = dof, EmbeddedNode = embeddedNode.Node, HostNode = embeddedNode.EmbeddedInElement.Nodes[i], Element = embeddedNode.EmbeddedInElement };
                        if (!superElementMap.ContainsKey(superElementDOF))
                        {
                            superElementMap.Add(superElementDOF, index);
                            index++;
                        }
                    }
                }

                //var independentEmbeddedDOFs = model.NodalDOFsDictionary[embeddedNode.Node.ID].Select(x => x.Key).Except(embeddedNode.DependentDOFs);
                //int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(embeddedNode.Node);

                //var independentEmbeddedDOFs = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement].Except(embeddedNode.DependentDOFs);

                foreach (var dof in independentEmbeddedDOFs)
                {
                    var superElementDOF = new SuperElementDOF() { DOF = dof, EmbeddedNode = embeddedNode.Node, HostNode = null, Element = null };
                    if (!superElementMap.ContainsKey(superElementDOF))
                    {
                        superElementMap.Add(superElementDOF, index);
                        index++;
                    }
                }
            }

            foreach (var node in embeddedElement.INodes.Except(e.EmbeddedNodes.Select(x => x.Node)))
            {
                int nodeOrderInEmbeddedElement = embeddedElement.INodes.IndexOf(node);
                var currentNodeDOFs = embeddedElement.IElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement];
                foreach (var dof in currentNodeDOFs)
                {
                    var superElementDOF = new SuperElementDOF() { DOF = dof, EmbeddedNode = node, HostNode = null, Element = null };
                    if (!superElementMap.ContainsKey(superElementDOF))
                    {
                        superElementMap.Add(superElementDOF, index);
                        index++;
                    }
                }
            }
        }

        private void CalculateTransformationMatrix()
        {
            var e = embeddedElement.IElementType as IEmbeddedElement;
            int row = 0;
            int col = 0;
            int totalRows = embeddedElement.IElementType.DOFEnumerator.GetDOFTypes(embeddedElement).SelectMany(x => x).Count();
            //var matrix = new double[totalRows, superElementMap.Count];
            var transformationMatrixOriginal = new Matrix2D(totalRows, superElementMap.Count);

            foreach (var embeddedNode in e.EmbeddedNodes)
            {
                var localTransformationMatrix = transformation.GetTransformationVector(embeddedNode);
                var localHostDOFs = transformation.GetDOFTypesOfHost(embeddedNode);
                int nodeOrderInEmbeddedElement = embeddedElement.INodes.IndexOf(embeddedNode.Node);
                var embeddedNodeDOFQuantity = embeddedElement.IElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement].Count;
                int dependentDOFs = transformation.GetDependentDOFTypes.Count;

                for (int i = 0; i < dependentDOFs; i++)
                {
                    col = 0;
                    for (int j = 0; j < localHostDOFs.Count; j++)
                    {
                        for (int k = 0; k < localHostDOFs[j].Count; k++)
                        {
                            var superelement = new SuperElementDOF() { DOF = localHostDOFs[j][k], Element = embeddedNode.EmbeddedInElement, EmbeddedNode = embeddedNode.Node, HostNode = embeddedNode.EmbeddedInElement.Nodes[j] };
                            transformationMatrixOriginal[40+row + i, superElementMap[superelement]] = localTransformationMatrix[i][col];
                            col++;
                        }
                    }
                }
                row += dependentDOFs;

                var independentEmbeddedDOFs = embeddedElement.IElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement].Except(embeddedNode.DependentDOFs).ToArray();
                for (int j = 0; j < independentEmbeddedDOFs.Length; j++)
                {
                    var superelement = new SuperElementDOF() { DOF = independentEmbeddedDOFs[j], Element = null, HostNode = null, EmbeddedNode = embeddedNode.Node };
                    transformationMatrixOriginal[40+row, superElementMap[superelement]] = 1;
                    row++;
                }
            }

            foreach (var node in embeddedElement.INodes.Except(e.EmbeddedNodes.Select(x => x.Node)))
            {
                int nodeOrderInEmbeddedElement = embeddedElement.INodes.IndexOf(node);
                var currentNodeDOFs = embeddedElement.IElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement];
                for (int j = 0; j < currentNodeDOFs.Count; j++)
                {
                    var superelement = new SuperElementDOF() { DOF = currentNodeDOFs[j], Element = null, HostNode = null, EmbeddedNode = node };
                    transformationMatrixOriginal[row-24, superElementMap[superelement]] = 1;
                    row++;
                }
            }

            //StreamWriter sw = File.CreateText(String.Format("TransformationMatrix{0}.txt", embeddedElement.ID));
            //for (int i = 0; i < totalRows; i++)
            //{
            //    var line = String.Empty;
            //    for (int j = 0; j < superElementMap.Count; j++)
            //        line += matrix[i,j].ToString() + ";";
            //    sw.WriteLine(line);
            //}
            //sw.Close();
            //transformationMatrixOriginal = new Matrix2D(matrix);
            transformationMatrix = transformationMatrixOriginal;
        }

        private void Initialize()
        {
            var e = embeddedElement.IElementType as IEmbeddedElement;
            if (e == null) return;
            if (e.EmbeddedNodes.Count == 0) return;

            InitializeMappings();
            CalculateTransformationMatrix();
        }

        //public IMatrix2D<double> GetTransformedMatrix(IMatrix2D<double> matrix)
        //{
        //    var e = embeddedElement.ElementType as IEmbeddedElement;
        //    //if (e == null || !isElementEmbedded) return matrix;
        //    if (e == null) return matrix;
        //    if (e.EmbeddedNodes.Count == 0) return matrix;

        //    return transformationMatrix.Transpose() * ((SymmetricMatrix2D<double>)matrix).ToMatrix2D() * transformationMatrix;
        //}

        public IMatrix2D GetTransformedMatrix(IMatrix2D matrix)
        {
            var e = embeddedElement.IElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return matrix;
            if (e.EmbeddedNodes.Count == 0) return matrix;

            return transformationMatrix.Transpose() * (Matrix2D)matrix * transformationMatrix;
            //return transformationMatrix.MultiplyTransposeThisTimesOtherTimesThis((Matrix2D)matrix);
        }

        public double[] GetTransformedDisplacementsVector(double[] vector)
        {
            var e = embeddedElement.IElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return vector;
            if (e.EmbeddedNodes.Count == 0) return vector;

            var result = new double[transformationMatrix.Rows];
            transformationMatrix.Multiply(new Vector(vector), result);
            return result;
        }

        public double[] GetTransformedForcesVector(double[] vector)
        {
            var e = embeddedElement.IElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return vector;
            if (e.EmbeddedNodes.Count == 0) return vector;

            Matrix2D transpose = transformationMatrix.Transpose();
            var result = new double[transpose.Rows];
            transpose.Multiply(new Vector(vector), result);
            return result;
        }

        public IList<IList<DOFType>> GetDOFTypes(IElement element)
        {
            //return element.ElementType.GetElementDOFTypes(element);

            var dofs = new List<IList<DOFType>>();
            INode currentNode = null;
            List<DOFType> nodeDOFs = null;

            foreach (var superElement in superElementMap)
            {
                INode node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;

                if (currentNode != node)
                {
                    if (nodeDOFs != null)
                        dofs.Add(nodeDOFs);
                    currentNode = node;
                    nodeDOFs = new List<DOFType>();
                }
                nodeDOFs.Add(superElement.Key.DOF);
            }
            if (nodeDOFs != null)
                dofs.Add(nodeDOFs);

            return dofs;
        }

        public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(IElement element)
        {
            //if (embeddedElement != element) throw new ArgumentException();

            var nodesDictionary = new Dictionary<INode, int>();
            int index = 0;
            foreach (var node in element.INodes)
            {
                nodesDictionary.Add(node, index);
                index++;
            }
            
            var dofs = new List<IList<DOFType>>();
            for (int i = 0; i < element.INodes.Count; i++)
                dofs.Add(new List<DOFType>());

            INode currentNode = null;
            List<DOFType> nodeDOFs = null;

            foreach (var superElement in superElementMap)
            {
                if (superElement.Key.HostNode != null) continue;
                INode node = superElement.Key.EmbeddedNode;
                //Node node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;

                if (currentNode != node)
                {
                    if (nodeDOFs != null)
                        dofs[nodesDictionary[currentNode]] = nodeDOFs;
                    currentNode = node;
                    nodeDOFs = new List<DOFType>();
                }
                nodeDOFs.Add(superElement.Key.DOF);
            }
            if (nodeDOFs != null)
                dofs[nodesDictionary[currentNode]] = nodeDOFs;
            //dofs.Add(nodeDOFs);

            return dofs;
        }

        public IList<INode> GetNodesForMatrixAssembly(IElement element)
        {
            var nodes = new List<INode>();
            INode currentNode = null;
            foreach (var superElement in superElementMap)
            {
                INode node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;
                if (currentNode != node)
                {
                    if (currentNode != null) 
                        nodes.Add(currentNode);
                    currentNode = node;
                }
                //if (nodes.IndexOf(node) < 0)
                //    nodes.Add(node);
            }
            if (currentNode != null)
                nodes.Add(currentNode);

            return nodes;
        }
    }
}
