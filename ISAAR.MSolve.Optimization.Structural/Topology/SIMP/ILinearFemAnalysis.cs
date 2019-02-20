using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: The methods depending on if the mesh is uniform, if it has passive elements etc should be in a different interface than
//      the methods concerning the execution of the analysis.
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP
{
    public interface ILinearFemAnalysis
    {
        int NumElements { get; }
        int NumLoadCases { get; }

        void AnalyzeModelWithDensities(IVectorView densities);

        // What about passive elements?
        double CalculateTotalVolume(IVectorView densities);

        Vector GetElementDisplacements(int elementIdx, int loadCaseIdx);

        /// <summary>
        /// The element's stiffness matrix after it has been scaled with respect to the penalized density.
        /// </summary>
        IMatrixView GetPenalizedStiffnessOfElement(int elementIdx, IVectorView densities);

        void Initialize();
    }
}
