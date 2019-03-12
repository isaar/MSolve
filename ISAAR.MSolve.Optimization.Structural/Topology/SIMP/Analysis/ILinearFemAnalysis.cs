using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: The methods depending on if the mesh is uniform, if it has passive elements etc should be in a different interface than
//      the methods concerning the execution of the analysis.
//TODO: Add an implementation for uniform meshes that caches the unit stiffness matrix. That class would not assume a hard coded
//      node and dof ordering. The assumption would be that elements have identical geometry, whatever that may be.
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.Analysis
{
    public interface ILinearFemAnalysis
    {
        int NumElements { get; }
        int NumLoadCases { get; }

        // What about passive elements?
        double CalculateTotalVolume(IVectorView densities);

        Vector GetElementDisplacements(int elementIdx, int loadCaseIdx);

        /// <summary>
        /// The element's stiffness matrix for the E that corresponds to it at the current topology step.
        /// </summary>
        IMatrixView GetElementStiffnessForCurrentMaterial(int elementIdx);

        /// <summary>
        /// The element's stiffness for E = 1.0
        /// </summary>
        IMatrixView GetElementStiffnessForUnitMaterial(int elementIdx);

        void Initialize();
        void UpdateModelAndAnalyze(IVectorView youngModuli);
    }
}
