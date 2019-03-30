using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Finite element for heat transfer along a single direction. Can be used for 1D, 2D and 3D problems. Does not take into 
    /// account geometric non-linearities.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ThermalRod : IFiniteElement_v2, IEmbeddedElement_v2
    {
        private const int numNodes = 2;
        private const int numDofs = 2;
        private static readonly DOFType[][] dofTypes = {
            new DOFType[] { DOFType.Temperature }, new DOFType[] { DOFType.Temperature } }; 

        private readonly ThermalMaterial material;

        public ThermalRod(IReadOnlyList<Node_v2> nodes, double crossSectionArea, ThermalMaterial material)
        {
            Debug.Assert(nodes.Count == 2, "Thermal rod element must have exactly 2 nodes.");
            this.material = material;
            this.Nodes = nodes;
            this.CrossSectionArea = crossSectionArea;
            this.Length = nodes[0].CalculateEuclidianDistanceFrom(nodes[1]);
        }

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public double CrossSectionArea { get; }
        public double Length { get; }
        public IReadOnlyList<Node_v2> Nodes { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDofEnumerator_v2 DofEnumerator { get; set; } = new GenericDofEnumerator_v2();

        public IList<EmbeddedNode_v2> EmbeddedNodes { get; } = new List<EmbeddedNode_v2>();

        public IMatrix MassMatrix(IElement_v2 element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix BuildCapacityMatrix()
        {
            double kdAL = material.SpecialHeatCoeff * material.Density * CrossSectionArea * Length;
            double[,] capacity = { { kdAL / 3.0, kdAL/ 6.0 }, { kdAL / 6.0, kdAL / 3.0 } };
            return Matrix.CreateFromArray(capacity);
        }

        public Matrix BuildConductivityMatrix()
        {

            double cAoverL = material.ThermalConductivity * CrossSectionArea / Length;
            double[,] conductivity = { { cAoverL, -cAoverL }, { -cAoverL, cAoverL } };
            return Matrix.CreateFromArray(conductivity);
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element) => dofTypes;

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element_v2 element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element_v2 element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element_v2 element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element_v2 element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        public IMatrix StiffnessMatrix(IElement_v2 element)
        {
            return DofEnumerator.GetTransformedMatrix(BuildConductivityMatrix());
        }

        public IMatrix DampingMatrix(IElement_v2 element)
        {
            throw new NotImplementedException();
        }

        public Dictionary<DOFType, int> GetInternalNodalDOFs(Element_v2 element, Node_v2 node)
        {
            if (node.ID == this.Nodes[0].ID) return new Dictionary<DOFType, int> { { DOFType.Temperature, 0 } };
            else if (node.ID == this.Nodes[1].ID) return new Dictionary<DOFType, int> { { DOFType.Temperature, 1 } };
            else throw new ArgumentException($"GetInternalNodalDOFs: Node {node.ID} not found in element {element.ID}.");
        }

        public double[] GetLocalDOFValues(Element_v2 hostElement, double[] hostDOFValues)
        {
            return DofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
        }
    }
}
