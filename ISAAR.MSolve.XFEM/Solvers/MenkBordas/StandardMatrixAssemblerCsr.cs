using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.XFEM.Solvers.MenkBordas
{
    class StandardMatrixAssemblerCsr : IStandardMatrixAssembler
    {
        public DokRowMajor Kss { get; set; }
        public DokRowMajor Ksc { get; set; }

        public void BuildStandardMatrices(Model2D model, XClusterDofOrderer globalDofOrderer)
        {
            int numDofsConstrained = globalDofOrderer.NumConstrainedDofs;
            int numDofsStandard = globalDofOrderer.NumStandardDofs;
            Kss = DokRowMajor.CreateEmpty(numDofsStandard, numDofsStandard);
            Ksc = DokRowMajor.CreateEmpty(numDofsStandard, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrix and add its contributions to the global matrices
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                globalDofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kss.AddSubmatrixSymmetric(kss, mapStandard);
                Ksc.AddSubmatrix(kss, mapStandard, mapConstrained);
            }
        }
    }
}
