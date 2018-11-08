//using ISAAR.MSolve.Analyzers.Interfaces;
//using ISAAR.MSolve.FEM.Entities;
//using ISAAR.MSolve.LinearAlgebra.Vectors;
//using System;
//using System.Collections.Generic;
//using System.Diagnostics;
//using System.Text;

//namespace ISAAR.MSolve.Analyzers
//{
//    public class NonLinearSubdomainUpdater_v2 : INonLinearSubdomainUpdater_v2
//    {
//        private readonly Subdomain subdomain;

//        public NonLinearSubdomainUpdater_v2(Subdomain subdomain)
//        {
//            this.subdomain = subdomain;
//        }

//        public IVector GetRHSFromSolution(IVectorView solution, IVectorView dSolution) //TODO leave 
//        {
//            return this.subdomain.GetRHSFromSolution_v2(solution, dSolution);
//        }

//        public void ResetState()
//        {
//            this.subdomain.ClearMaterialStresses();
//        }

//        public void UpdateState()
//        {
//            this.subdomain.SaveMaterialState();
//        }
//    }
//}

