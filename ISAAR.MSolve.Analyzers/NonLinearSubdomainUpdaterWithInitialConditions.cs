using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace ISAAR.MSolve.Analyzers
{
    //public class NonLinearSubdomainUpdaterWithInitialConditions : INonLinearSubdomainUpdater
    //{
    //    private readonly Subdomain subdomain;

    //    public NonLinearSubdomainUpdaterWithInitialConditions(Subdomain subdomain)
    //    {
    //        this.subdomain = subdomain;
    //    }

    //    public IVector GetRHSFromSolutionWithInitialDisplacemntsEffect(IVector solution, IVector dSolution, Dictionary<int, Node> boundaryNodes,
    //        Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements, Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements,
    //        int nIncrement, int totalIncrements) //TODO leave 
    //    {
    //        return this.subdomain.GetRHSFromSolutionWithInitialDisplacemntsEffect( solution, dSolution, boundaryNodes,
    //         initialConvergedBoundaryDisplacements,  totalBoundaryDisplacements,
    //         nIncrement, totalIncrements);
    //    }

    //    public void ResetState()
    //    {
    //        this.subdomain.ClearMaterialStresses();
    //    }

    //    public void UpdateState()
    //    {
    //        this.subdomain.SaveMaterialState();
    //    }

    //    public void ScaleConstraints(double scalingFactor)
    //    {
    //        throw new NotSupportedException();
    //    }

    //    public IVector GetRHSFromSolution(IVector solution, IVector dSolution) //TODO leave 
    //    {
    //        throw new NotSupportedException();
    //        return this.subdomain.GetRHSFromSolution(solution, dSolution);
    //    }
    //}
}
