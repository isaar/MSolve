using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Numerical.Commons;

//TODO: remove code that calculates rhs vector components (nodal loads, constraints, etc). It should be moved to dedicated 
//      classes like EquivalentLoadAssembler, so that it can be reused between subdomains of different projects (FEM, IGA, XFEM).
//TODO: same for multiscale
namespace ISAAR.MSolve.FEM.Entities
{
    public class Subdomain_v2 : ISubdomain_v2
    {
        private readonly List<Node_v2> nodes = new List<Node_v2>();

        public Subdomain_v2(int id)
        {
            this.ID = id;
        }

        public Table<INode, DOFType, double> Constraints { get; } = new Table<INode, DOFType, double>();

        IReadOnlyList<IElement_v2> ISubdomain_v2.Elements => Elements;
        public List<Element_v2> Elements { get; } = new List<Element_v2>();

        //public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

        public int ID { get; }

        IReadOnlyList<INode> ISubdomain_v2.Nodes => nodes;
        public IReadOnlyList<Node_v2> Nodes => nodes;

        public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }
        public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

        public Vector Forces { get; set; } //TODO: this doesn't belong here

        //public bool MaterialsModified
        //{
        //    get
        //    {
        //        bool modified = false;
        //        foreach (Element element in elementsDictionary.Values)
        //            if (element.ElementType.MaterialModified)
        //            {
        //                modified = true;
        //                break;
        //            }
        //        return modified;
        //    }
        //}
        public bool MaterialsModified { get; set; }
        
        //TODO: This belongs in EquivalentLoadsAssembler
        //TODO: the constraintScalingFactor parameter is not used.
        public double[] CalculateElementIncrementalConstraintDisplacements(IElement_v2 element, double constraintScalingFactor)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
            ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public double[] CalculateElementDisplacements(Element_v2 element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public void ClearMaterialStresses()
        {
            foreach (Element_v2 element in Elements) element.ElementType.ClearMaterialStresses();
        }

        public void DefineNodesFromElements()
        {
            nodes.Clear();
            var nodeComparer = Comparer<Node_v2>.Create((Node_v2 node1, Node_v2 node2) => node1.ID - node2.ID);
            var nodeSet = new SortedSet<Node_v2>(nodeComparer);
            foreach (Element_v2 element in Elements)
            {
                foreach (Node_v2 node in element.Nodes) nodeSet.Add(node);
            }
            nodes.AddRange(nodeSet);

            //foreach (var e in modelEmbeddedNodes.Where(x => nodeIDs.IndexOf(x.Node.ID) >= 0))
            //    EmbeddedNodes.Add(e);
        }

        //TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
        //      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
        //      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
        //      It is too easy to access the wrong instance of the constraint. 
        public void ExtractConstraintsFromGlobal(Table<INode, DOFType, double> globalConstraints)
        {
            //TODO: perhaps it is more efficient to traverse the global constraints instead of the subdomain's nodes, provided
            //      the latter are stored as a set. 
            //TODO: the next could be a Table method: Table.KeepDataOfRows(IEnumerable<TRow> rows)
            foreach (Node_v2 node in Nodes)
            {
                bool isNodeConstrained = globalConstraints.TryGetDataOfRow(node,
                    out IReadOnlyDictionary<DOFType, double> constraintsOfNode);
                if (isNodeConstrained)
                {
                    foreach (var dofDisplacementPair in constraintsOfNode)
                    {
                        Constraints[node, dofDisplacementPair.Key] = dofDisplacementPair.Value;
                    }
                }
            }

            // This is probably faster but assumes that nodes store their prescribed displacements, which I hate.
            //foreach (Node node in Nodes)
            //{
            //    if (node.Constraints == null) continue;
            //    foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
            //}
        }

        public int[] GetCornerNodes()
        {
            int nodex1y1z1 = -1;
            int nodex2y1z1 = -1;
            int nodex1y2z1 = -1;
            int nodex2y2z1 = -1;
            int nodex1y1z2 = -1;
            int nodex2y1z2 = -1;
            int nodex1y2z2 = -1;
            int nodex2y2z2 = -1;
            double x1 = Double.MaxValue;
            double x2 = Double.MinValue;
            double y1 = Double.MaxValue;
            double y2 = Double.MinValue;
            double z1 = Double.MaxValue;
            double z2 = Double.MinValue;

            foreach (var kv in nodes)
            {
                if (x1 > kv.X) x1 = kv.X;
                if (x2 < kv.X) x2 = kv.X;
                if (y1 > kv.Y) y1 = kv.Y;
                if (y2 < kv.Y) y2 = kv.Y;
                if (z1 > kv.Z) z1 = kv.Z;
                if (z2 < kv.Z) z2 = kv.Z;

                if (x1 == kv.X && y1 == kv.Y && z1 == kv.Z) nodex1y1z1 = kv.ID;
                if (x2 == kv.X && y1 == kv.Y && z1 == kv.Z) nodex2y1z1 = kv.ID;
                if (x1 == kv.X && y2 == kv.Y && z1 == kv.Z) nodex1y2z1 = kv.ID;
                if (x2 == kv.X && y2 == kv.Y && z1 == kv.Z) nodex2y2z1 = kv.ID;
                if (x1 == kv.X && y1 == kv.Y && z2 == kv.Z) nodex1y1z2 = kv.ID;
                if (x2 == kv.X && y1 == kv.Y && z2 == kv.Z) nodex2y1z2 = kv.ID;
                if (x1 == kv.X && y2 == kv.Y && z2 == kv.Z) nodex1y2z2 = kv.ID;
                if (x2 == kv.X && y2 == kv.Y && z2 == kv.Z) nodex2y2z2 = kv.ID;
            }

            return new[] { nodex1y1z1, nodex2y1z1, nodex1y2z1, nodex2y2z1, nodex1y1z2, nodex2y1z2, nodex1y2z2, nodex2y2z2 };
        }

        public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
        {
            var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element_v2 element in Elements)
            {
                //var localSolution = GetLocalVectorFromGlobal(element, solution);//TODOMaria: This is where the element displacements are calculated //removeMaria
                //var localdSolution = GetLocalVectorFromGlobal(element, dSolution);//removeMaria

                //TODO: ElementType should operate with Vector instead of double[]. Then the ToRawArray() calls can be removed
                double[] localSolution = CalculateElementDisplacements(element, solution);
                double[] localdSolution = CalculateElementDisplacements(element, dSolution);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain.MaterialsModified = true;
                var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public void ResetMaterialsModifiedProperty()
        {
            this.MaterialsModified = false;
            foreach (Element_v2 element in Elements) element.ElementType.ResetMaterialModified();
        }

        public void SaveMaterialState()
        {
            foreach (Element_v2 element in Elements) element.ElementType.SaveMaterialState();
        }

        //TODO: I am against modifying the constraints table of the subdomain. Instead the analyzer should keep a constraint
        //      displacements vector at global/subdomain scale and modify that.
        public void ScaleConstraints(double scalingFactor) => Constraints.ModifyValues((u) => scalingFactor * u);

        //TODO: this should probably belong to a dedicated class, which is abstracted by an interface. That class can be reused
        //      by various analyzer classes and as such does not need rewriting for IGA, XFEM or other problem types. Also it 
        //      can perform optimizations, such as using a prescribed displacements vector and element constrained dof maps,
        //      similarly to how free displacements are handled by ISubdomainFreeDofOrdering
        private static void ApplyConstraintDisplacements(IElement_v2 element, double[] elementNodalDisplacements,
            Table<INode, DOFType, double> constraints)
        {
            int elementDofIdx = 0;
            IList<INode> nodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IList<IList<DOFType>> dofs = element.ElementType.DofEnumerator.GetDOFTypes(element);
            for (int i = 0; i < nodes.Count; ++i)
            {
                //bool isConstrainedNode = constraintsDictionary.TryGetValue(nodes[i].ID, 
                //    out Dictionary<DOFType, double> constrainedDOFs);
                bool isConstrainedNode = constraints.TryGetDataOfRow(nodes[i],
                    out IReadOnlyDictionary<DOFType, double> constrainedDOFs);
                if (isConstrainedNode)
                {
                    foreach (DOFType dofType in dofs[i])
                    {
                        bool isConstrainedDof = constrainedDOFs.TryGetValue(dofType, out double constraintDisplacement);
                        //if (isConstrainedNode && isConstrainedDof)
                        if (isConstrainedDof)
                        {
                            Debug.Assert(elementNodalDisplacements[elementDofIdx] == 0); // TODO: and why is this an assumption?
                            elementNodalDisplacements[elementDofIdx] = constraintDisplacement;
                        }
                        ++elementDofIdx;
                    }
                }
                else elementDofIdx += dofs[i].Count;
            }
        }

        //ADDED1
        public IVector GetRHSFromSolutionWithInitialDisplacemntsEffect(IVectorView solution, IVectorView dSolution, Dictionary<int, Node_v2> boundaryNodes,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {
            ////prosthiki print

            //ekteleseis_counter2 += 1;
            //string counter_data = ekteleseis_counter2.ToString();
            //string path = string.Format(string2, counter_data);
            ////solution.WriteToFile(path);
            //double[] solution_data = new double[solution.Length];
            //solution_data = solution.CopyToArray();
            //WriteToFileVector(solution_data, path);

            var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element_v2 element in Elements)
            {
                var localSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, solution);
                ImposePrescribedDisplacementsWithInitialConditionSEffect(element, localSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
                var localdSolution = GetLocalVectorFromGlobalWithoutPrescribedDisplacements(element, dSolution);
                ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(element, localdSolution, boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, nIncrement, totalIncrements);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Subdomain.MaterialsModified = true;
                var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public double[] GetLocalVectorFromGlobalWithoutPrescribedDisplacements(Element_v2 element, IVectorView globalDisplacementVector)
        {
            double[] elementNodalDisplacements = FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            return elementNodalDisplacements;
        }
        

        public void ImposePrescribedDisplacementsWithInitialConditionSEffect(Element_v2 element, double[] localSolution, Dictionary<int, Node_v2> boundaryNodes,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {

            var elementDOFTypes = element.ElementType.DofEnumerator.GetDOFTypes(element);
            var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            int iElementMatrixColumn = 0;
            for (int j = 0; j < elementDOFTypes.Count; j++)
            {
                INode nodeColumn = matrixAssemblyNodes[j];
                int nodalDofsNumber = elementDOFTypes[j].Count;
                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                {
                    Dictionary<DOFType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
                    Dictionary<DOFType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
                    int positionOfDofInNode = 0;
                    foreach (DOFType doftype1 in elementDOFTypes[j])
                    {
                        if (nodalConvergedDisplacements.ContainsKey(doftype1))
                        {
                            localSolution[iElementMatrixColumn + positionOfDofInNode] = nodalConvergedDisplacements[doftype1] + (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
                            // TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
                        }
                        positionOfDofInNode += 1;
                    }
                }
                iElementMatrixColumn += nodalDofsNumber;
            }

        }

        public void ImposePrescribed_d_DisplacementsWithInitialConditionSEffect(Element_v2 element, double[] localSolution, Dictionary<int, Node_v2> boundaryNodes,
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
            int nIncrement, int totalIncrements)
        {

            var elementDOFTypes = element.ElementType.DofEnumerator.GetDOFTypes(element);
            var matrixAssemblyNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            int iElementMatrixColumn = 0;
            for (int j = 0; j < elementDOFTypes.Count; j++)
            {
                INode nodeColumn = matrixAssemblyNodes[j];
                int nodalDofsNumber = elementDOFTypes[j].Count;
                if (boundaryNodes.ContainsKey(nodeColumn.ID))
                {
                    Dictionary<DOFType, double> nodalConvergedDisplacements = initialConvergedBoundaryDisplacements[nodeColumn.ID];
                    Dictionary<DOFType, double> nodalTotalDisplacements = totalBoundaryDisplacements[nodeColumn.ID];
                    int positionOfDofInNode = 0;
                    foreach (DOFType doftype1 in elementDOFTypes[j])
                    {
                        if (nodalConvergedDisplacements.ContainsKey(doftype1))
                        {
                            localSolution[iElementMatrixColumn + positionOfDofInNode] = (nodalTotalDisplacements[doftype1] - nodalConvergedDisplacements[doftype1]) * ((double)nIncrement / (double)totalIncrements);
                            // 1) den vazoume mono (1/increments) alla (nIncrement/increments) dioti metaxu aftwn twn nIncrements den exei mesolavhsei save sta material ths mikroklimakas
                            // TODO: this can be done faster: create a dictionary<...,dictionary> with the difference of the two values and use that and precalculate coefficient for scaling
                        }
                        positionOfDofInNode += 1;
                    }
                }
                iElementMatrixColumn += nodalDofsNumber;
            }

        }
    }
}
