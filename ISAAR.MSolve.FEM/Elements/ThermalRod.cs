using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Finite element for heat transfer along a single direction. Can be used for 1D, 2D and 3D problems. Does not take into 
    /// account geometric non-linearities.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ThermalRod : IFiniteElement
    {
        private const int numNodes = 2;
        private const int numDofs = 2;
        private static readonly DOFType[][] dofTypes = {
            new DOFType[] { DOFType.Temperature }, new DOFType[] { DOFType.Temperature } }; 

        private readonly ThermalMaterial material;

        public ThermalRod(IReadOnlyList<Node2D> nodes, double crossSectionArea, ThermalMaterial material)
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
        public IReadOnlyList<Node2D> Nodes { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDOFEnumerator DOFEnumerator { get; set; } = new GenericDOFEnumerator();

        public IMatrix2D MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix2D BuildCapacityMatrix()
        {
            double kdAL = material.SpecialHeatCoeff * material.Density * CrossSectionArea * Length;
            double[,] capacity = { { kdAL / 3.0, kdAL/ 6.0 }, { kdAL / 6.0, kdAL / 3.0 } };
            return new Matrix2D(capacity);
        }

        public Matrix2D BuildConductivityMatrix()
        {

            double cAoverL = material.ThermalConductivity * CrossSectionArea / Length;
            double[,] conductivity = { { cAoverL, -cAoverL }, { -cAoverL, cAoverL } };
            return new Matrix2D(conductivity);
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element) => dofTypes;

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
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

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            return BuildConductivityMatrix();
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }
    }
}
