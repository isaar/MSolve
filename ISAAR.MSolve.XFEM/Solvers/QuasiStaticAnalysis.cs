using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.CrackGeometry;
using ISAAR.MSolve.XFEM.CrackPropagation;
using ISAAR.MSolve.XFEM.CrackPropagation.Direction;
using ISAAR.MSolve.XFEM.CrackPropagation.Jintegral;
using ISAAR.MSolve.XFEM.CrackPropagation.Length;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Mesh;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class QuasiStaticAnalysis
    {
        private readonly Model2D model;
        private readonly IMesh2D<XNode2D, XContinuumElement2D> mesh;
        private readonly IExteriorCrack crack;
        private readonly ISolver solver;
        private readonly Propagator propagator;
        private readonly double fractureToughness;
        private readonly int maxIterations;

        public QuasiStaticAnalysis(Model2D model, IMesh2D<XNode2D, XContinuumElement2D> mesh, IExteriorCrack crack,
            ISolver solver, Propagator propagator, double fractureToughness, int maxIterations)
        {
            this.model = model;
            this.mesh = mesh;
            this.crack = crack;
            this.solver = solver;
            this.propagator = propagator;
            this.fractureToughness = fractureToughness;
            this.maxIterations = maxIterations;
        }

        public IReadOnlyList<ICartesianPoint2D> Analyze()
        {
            var crackPath = new List<ICartesianPoint2D>();
            crackPath.Add(crack.CrackMouth);
            crackPath.Add(crack.GetCrackTip(CrackTipPosition.Single));

            int iteration;
            try
            {
                solver.Initialize();
                for (iteration = 0; iteration < maxIterations; ++iteration)
                {
                    //Console.WriteLine(
                        //"********************************** Iteration {0} **********************************", iteration);
                    crack.UpdateEnrichments();
                    solver.Solve();

                    Vector totalConstrainedDisplacements = model.CalculateConstrainedDisplacements(solver.DOFEnumerator);
                    Vector totalFreeDisplacements = solver.Solution;

                    double growthAngle, growthIncrement;
                    propagator.Propagate(solver.DOFEnumerator, totalFreeDisplacements, totalConstrainedDisplacements,
                        out growthAngle, out growthIncrement);
                    crack.UpdateGeometry(growthAngle, growthIncrement);
                    ICartesianPoint2D newTip = crack.GetCrackTip(CrackTipPosition.Single);
                    crackPath.Add(newTip);

                    double sifEffective = EquivalentSIF(propagator.Logger.SIFsMode1[iteration],
                        propagator.Logger.SIFsMode2[iteration]);
                    if (sifEffective >= fractureToughness)
                    {
                        Console.WriteLine(
                            "Propagation analysis terminated: Failure due to fracture tougness being exceeded.");
                        break;
                    }
                    else if (!mesh.IsInsideBoundary(newTip))
                    {
                        Console.WriteLine(
                            "Propagation analysis terminated: Failure due to the crack reaching the domain's boudary.");
                        break;
                    }
                }
                if (iteration == maxIterations)
                {
                    Console.WriteLine(
                    "Propagation analysis terminated: All {0} iterations were completed.", maxIterations);
                }
            }
            catch(Exception ex)
            {
                Console.WriteLine("Crack path: X Y");
                foreach (var point in crackPath)
                {
                    Console.WriteLine("{0} {1}", point.X, point.Y);
                }
                throw ex;
            }

            return crackPath;
        }

        // TODO: Abstract this and add Tanaka_1974 approach
        private double EquivalentSIF(double sifMode1, double sifMode2)
        {
            return Math.Sqrt(sifMode1 * sifMode1 + sifMode2 * sifMode2);
        }

    }
}
