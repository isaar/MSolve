using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti1;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Feti.Feti1
{
    /// <summary>
    /// Tests from Papagiannakis bachelor thesis (NTUA 2011), p. 97 - 100
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class HomogeneousCantileverTests
    {
        
    }
}
