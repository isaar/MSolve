using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    public class NonLinearSubdomainUpdater : INonLinearSubdomainUpdater
    {
        private readonly Subdomain subdomain;

        public NonLinearSubdomainUpdater(Subdomain subdomain)
        {
            this.subdomain = subdomain;
        }

        public IVector GetRHSFromSolution(IVector solution, IVector dSolution) //TODO leave 
        {
            return this.subdomain.GetRHSFromSolution(solution, dSolution);
        }

        public void ResetState()
        {
            this.subdomain.ClearMaterialStresses();
        }

        public void UpdateState()
        {
            this.subdomain.SaveMaterialState();
        }
    }
}

