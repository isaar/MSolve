using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using Moq;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Tests
{
    public class Beam2DNewmarkDynamicAanalysisTest
    {
        [Fact]
        public void LinearElasticBeam2DNewmarkDynamicAnalysisTest()
        {
            double youngModulus = 2.0e08;
            double poissonRatio = 0.3;
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });
            int totalElements = 10;
            int totalNodes = totalElements + 1;
            ElasticMaterial material = new ElasticMaterial()
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };
            int iNodeID = 0;
            for (int iNode = 0; iNode < totalNodes; iNode++)
            {
                iNodeID++;
                var n = new Node() { ID = iNodeID, X = iNodeID * 10, Y = 0, Z = 0 };
                model.NodesDictionary.Add(iNodeID, n);
            }
            iNodeID = 0;
            int iElementID = 0;
            for (int iElement = 0; iElement < totalElements; iElement++)
            {
                iElementID++;
                iNodeID++;
                double b = 0.30;
                double h = 0.60;
                var e = new Element()
                {
                    ID = iElementID,
                    ElementType = new EulerBeam2D(youngModulus)
                    {
                        Density = 7.85,
                        SectionArea = b * h,
                        MomentOfInertia = (1/12) * b * h * h * h,
                    }
                };
                e.NodesDictionary.Add(iNodeID, model.NodesDictionary[iNodeID]);
                e.NodesDictionary.Add(iNodeID + 1, model.NodesDictionary[iNodeID + 1]);                
                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[1].ElementsDictionary.Add(e.ID, e);
            }
            model.Loads.Add(new Load() { Amount = -100, Node = model.NodesDictionary[totalNodes], DOF = DOFType.Y });
            model.ConnectDataStructures();
            var linearSystems = new Dictionary<int, ILinearSystem>();
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            SolverSkyline solver = new SolverSkyline(linearSystems[1]);
            ProblemStructural provider = new ProblemStructural(model, linearSystems);
            LinearAnalyzer analyzer = new LinearAnalyzer(solver, linearSystems);
            NewmarkDynamicAnalyzer dynamicAnalyzer = new NewmarkDynamicAnalyzer(provider, analyzer, linearSystems, 0.25, 0.5, 0.28, 3.36);
            dynamicAnalyzer.BuildMatrices();
            dynamicAnalyzer.Initialize();
            dynamicAnalyzer.Solve();
            Assert.Equal(2.2840249264795207, linearSystems[1].Solution[0], 8);
            Assert.Equal(2.4351921891904156, linearSystems[1].Solution[1], 8);
        }
    }
}
