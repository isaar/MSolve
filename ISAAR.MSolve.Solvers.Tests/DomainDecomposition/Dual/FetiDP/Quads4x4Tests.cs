using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.InterfaceProblem;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Pcg;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Preconditioning;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.StiffnessDistribution;
using Xunit;

//TODO: Also test stiffness distribution and preconditioners in other classes.
//TODO: Create the dofSeparator and lagrangeEnumerator manually, without using FetiDPSolver.
//TODO: TestInterfaceProblemSolution should mock matrices and vectors from TestInterfaceProblemCreation.
//TODO: There is a lot of code duplication between the methods.
namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Quads4x4Tests
    {
        private static Dictionary<int, Matrix> MatricesKrr
        {
            get
            {
                var Krr = new Dictionary<int, Matrix>();
                Krr[0] = Matrix.CreateFromArray(new double[,]
                {
                    { 0.9890109890, 0, 0.10989010990, 0, -0.24725274730, -0.17857142860, 0, 0 },
                    { 0, 0.9890109890,  0, -0.60439560440, -0.17857142860, -0.24725274730, 0, 0 },
                    { 0.10989010990, 0, 1.9780219780,  0, -0.60439560440, 0, 0.10989010990, 0 },
                    { 0, -0.60439560440, 0, 1.9780219780,  0, 0.10989010990, 0, -0.60439560440 },
                    { -0.24725274730,  -0.17857142860, -0.60439560440, 0, 0.9890109890,  0, -0.24725274730, 0.17857142860 },
                    { -0.17857142860,  -0.24725274730, 0, 0.10989010990, 0, 0.9890109890,  0.17857142860, -0.24725274730 },
                    { 0, 0, 0.10989010990, 0, -0.24725274730, 0.17857142860, 0.9890109890,  0 },
                    { 0, 0, 0, -0.60439560440, 0.17857142860, -0.24725274730, 0, 0.9890109890 }
                });

                Krr[1] = Matrix.CreateFromArray(new double[,]
                {
                    { 0.9890109890,    0,  -0.3021978022, -0.01373626370,-0.2472527473, 0.1785714286,  0.1098901099,  0,  -0.2472527473, -0.1785714286, 0,  0      },
                    { 0, 0.9890109890,  0.01373626370, 0.05494505490, 0.1785714286,  -0.2472527473, 0,  -0.6043956044, -0.1785714286, -0.2472527473, 0,  0         },
                    { -0.3021978022,   0.01373626370, 0.4945054945,  -0.1785714286, 0,  0,  -0.2472527473, 0.1785714286,  0.05494505490, -0.01373626370,    0,  0  },
                    { -0.01373626370,  0.05494505490, -0.1785714286, 0.4945054945,  0,  0,  0.1785714286,  -0.2472527473, 0.01373626370, -0.3021978022, 0,  0      },
                    { -0.2472527473,   0.1785714286,  0,  0,  0.9890109890,  0,  -0.6043956044, 0,  0,  0,  -0.2472527473, -0.1785714286                           },
                    { 0.1785714286,    -0.2472527473, 0,  0,  0,  0.9890109890,  0,  0.1098901099,  0,  0,  -0.1785714286, -0.2472527473                           },
                    { 0.1098901099,    0,  -0.2472527473, 0.1785714286,  -0.6043956044, 0,  1.978021978,   0,  -0.6043956044, 0,  0.1098901099,  0                 },
                    { 0, -0.6043956044, 0.1785714286,  -0.2472527473, 0,  0.1098901099,  0,  1.978021978,   0,  0.1098901099,  0,  -0.6043956044                   },
                    { -0.2472527473,   -0.1785714286, 0.05494505490, 0.01373626370, 0,  0,  -0.6043956044, 0,  0.9890109890,  0,  -0.2472527473, 0.1785714286      },
                    { -0.1785714286,   -0.2472527473, -0.01373626370,    -0.3021978022, 0,  0,  0,  0.1098901099,  0,  0.9890109890,  0.1785714286,  -0.2472527473 },
                    { 0,    0,  0,  0,  -0.2472527473, -0.1785714286, 0.1098901099,  0,  -0.2472527473, 0.1785714286,  0.9890109890,  0                            },
                    { 0,    0,  0,  0,  -0.1785714286, -0.2472527473, 0,  -0.6043956044, 0.1785714286,  -0.2472527473, 0,  0.9890109890                            }
                });

                Krr[2] = Matrix.CreateFromArray(new double[,]
                {
                    { 0.9890109890, 0, 0.10989010990, 0, -0.24725274730, -0.17857142860, 0, 0 },
                    { 0, 0.9890109890, 0, -0.60439560440, -0.17857142860, -0.24725274730, 0, 0 },
                    { 0.10989010990, 0, 1.9780219780, 0, -0.60439560440, 0, 0.10989010990, 0 },
                    { 0, -0.60439560440, 0, 1.9780219780, 0, 0.10989010990, 0, -0.60439560440 },
                    { -0.24725274730, -0.17857142860, -0.60439560440, 0, 0.9890109890, 0, -0.24725274730, 0.17857142860 },
                    { -0.17857142860, -0.24725274730, 0, 0.10989010990, 0, 0.9890109890, 0.17857142860, -0.24725274730 },
                    { 0, 0, 0.10989010990, 0, -0.24725274730, 0.17857142860, 0.9890109890, 0 },
                    { 0, 0, 0, -0.60439560440, 0.17857142860, -0.24725274730, 0, 0.9890109890 }
                });

                Krr[3] = Matrix.CreateFromArray(new double[,]
                {
                    { 0.9890109890, 0, -0.24725274730, 0.17857142860, 0.10989010990, 0, -0.24725274730, -0.17857142860, 0, 0, 0, 0                           },
                    { 0, 0.9890109890, 0.17857142860, -0.24725274730, 0, -0.60439560440, -0.17857142860, -0.24725274730, 0, 0, 0, 0                          },
                    { -0.24725274730, 0.17857142860, 0.9890109890, 0, -0.60439560440, 0, 0, 0, -0.24725274730, -0.17857142860, 0, 0                          },
                    { 0.17857142860, -0.24725274730, 0, 0.9890109890, 0, 0.10989010990, 0, 0, -0.17857142860, -0.24725274730, 0, 0                           },
                    { 0.10989010990, 0, -0.60439560440, 0, 1.9780219780, 0, -0.60439560440, 0, 0.10989010990, 0, -0.24725274730, -0.17857142860              },
                    { 0, -0.60439560440, 0, 0.10989010990, 0, 1.9780219780, 0, 0.10989010990, 0, -0.60439560440, -0.17857142860, -0.24725274730              },
                    { -0.24725274730, -0.17857142860, 0, 0, -0.60439560440, 0, 0.9890109890, 0, -0.24725274730, 0.17857142860, 0.05494505490, -0.01373626370 },
                    { -0.17857142860, -0.24725274730, 0, 0, 0, 0.10989010990, 0, 0.9890109890, 0.17857142860, -0.24725274730, 0.01373626370, -0.30219780220  },
                    { 0, 0, -0.24725274730, -0.17857142860, 0.10989010990, 0, -0.24725274730, 0.17857142860, 0.9890109890, 0, -0.30219780220, 0.01373626370  },
                    { 0, 0, -0.17857142860, -0.24725274730, 0, -0.60439560440, 0.17857142860, -0.24725274730, 0, 0.9890109890, -0.01373626370, 0.05494505490 },
                    { 0, 0, 0, 0, -0.24725274730, -0.17857142860, 0.05494505490, 0.01373626370, -0.30219780220, -0.01373626370, 0.49450549450, 0.17857142860 },
                    { 0, 0, 0, 0, -0.17857142860, -0.24725274730, -0.01373626370, -0.30219780220, 0.01373626370, 0.05494505490, 0.17857142860, 0.49450549450 }
                });
                return Krr;
            }
        }

        private static Dictionary<int, Matrix> MatricesKrc
        {
            get
            {
                var Krc = new Dictionary<int, Matrix>();
                Krc[0] = Matrix.CreateFromArray(new double[,]
                {
                    {-0.30219780220, -0.01373626370, 0, 0                         },
                    {0.01373626370, 0.05494505490, 0, 0                           },
                    {-0.24725274730, 0.17857142860, -0.24725274730, -0.17857142860},
                    {0.17857142860, -0.24725274730, -0.17857142860, -0.24725274730},
                    {0.05494505490, 0.01373626370, 0.05494505490, -0.01373626370  },
                    {-0.01373626370, -0.30219780220, 0.01373626370, -0.30219780220},
                    {0, 0, -0.30219780220, 0.01373626370                          },
                    {0, 0, -0.01373626370, 0.05494505490                          }
                });

                Krc[1] = Matrix.CreateFromArray(new double[,]
                {
                    {-0.30219780220, 0.01373626370, 0, 0, 0, 0                                                     },
                    {-0.01373626370, 0.05494505490, 0, 0, 0, 0                                                     },
                    {0, 0, 0, 0, 0, 0                                                                              },
                    {0, 0, 0, 0, 0, 0                                                                              },
                    {0.05494505490, -0.01373626370, 0.05494505490, 0.01373626370, 0, 0                             },
                    {0.01373626370, -0.30219780220, -0.01373626370, -0.30219780220, 0, 0                           },
                    {-0.24725274730, -0.17857142860, -0.24725274730, 0.17857142860, -0.24725274730, -0.17857142860 },
                    {-0.17857142860, -0.24725274730, 0.17857142860, -0.24725274730, -0.17857142860, -0.24725274730 },
                    {0, 0, 0, 0, 0.05494505490, -0.01373626370                                                     },
                    {0, 0, 0, 0, 0.01373626370, -0.30219780220                                                     },
                    {0, 0, -0.30219780220, -0.01373626370, -0.30219780220, 0.01373626370                           },
                    {0, 0, 0.01373626370, 0.05494505490, -0.01373626370, 0.05494505490                             }
                });

                Krc[2] = Matrix.CreateFromArray(new double[,]
                {
                    {-0.30219780220, -0.01373626370, 0, 0                          },
                    {0.01373626370, 0.05494505490, 0, 0                            },
                    {-0.24725274730, 0.17857142860, -0.24725274730, -0.17857142860 },
                    {0.17857142860, -0.24725274730, -0.17857142860, -0.24725274730 },
                    {0.05494505490, 0.01373626370, 0.05494505490, -0.01373626370   },
                    {-0.01373626370, -0.30219780220, 0.01373626370, -0.30219780220 },
                    {0, 0, -0.30219780220, 0.01373626370                           },
                    {0, 0, -0.01373626370, 0.05494505490                           }
                });

                Krc[3] = Matrix.CreateFromArray(new double[,]
                {
                    {-0.30219780220, 0.01373626370, -0.30219780220, -0.01373626370, 0, 0                         },
                    {-0.01373626370, 0.05494505490, 0.01373626370, 0.05494505490, 0, 0                           },
                    {0.05494505490, -0.01373626370, 0, 0, 0.05494505490, 0.01373626370                           },
                    {0.01373626370, -0.30219780220, 0, 0, -0.01373626370, -0.30219780220                         },
                    {-0.24725274730, -0.17857142860, -0.24725274730, 0.17857142860, -0.24725274730, 0.17857142860},
                    {-0.17857142860, -0.24725274730, 0.17857142860, -0.24725274730, 0.17857142860, -0.24725274730},
                    {0, 0, 0.05494505490, 0.01373626370, 0, 0                                                    },
                    {0, 0, -0.01373626370, -0.30219780220, 0, 0                                                  },
                    {0, 0, 0, 0, -0.30219780220, -0.01373626370                                                  },
                    {0, 0, 0, 0, 0.01373626370, 0.05494505490                                                    },
                    {0, 0, 0, 0, 0, 0                                                                            },
                    {0, 0, 0, 0, 0, 0                                                                            }
                });
                return Krc;
            }
        }

        private static Dictionary<int, Matrix> MatricesKcc
        {
            get
            {
                var Kcc = new Dictionary<int, Matrix>();
                Kcc[0] = Matrix.CreateFromArray(new double[,]
                {
                    {0.494505494500000, -0.178571428600000, 0, 0},
                    {-0.178571428600000, 0.494505494500000, 0, 0},
                    {0, 0, 0.494505494500000, 0.178571428600000},
                    {0, 0, 0.178571428600000, 0.494505494500000}
                });

                Kcc[1] = Matrix.CreateFromArray(new double[,]
                {
                    {0.494505494500000, 0.178571428600000, 0, 0, 0, 0},
                    {0.178571428600000, 0.494505494500000, 0, 0, 0, 0},
                    {0, 0, 0.494505494500000, -0.178571428600000, 0, 0},
                    {0, 0, -0.178571428600000, 0.494505494500000, 0, 0},
                    {0, 0, 0, 0, 0.494505494500000, 0.178571428600000},
                    {0, 0, 0, 0, 0.178571428600000, 0.494505494500000}
                });

                Kcc[2] = Matrix.CreateFromArray(new double[,]
                {
                    {0.494505494500000, -0.178571428600000, 0, 0},
                    {-0.178571428600000, 0.494505494500000, 0, 0},
                    {0, 0, 0.494505494500000, 0.178571428600000 },
                    {0, 0, 0.178571428600000, 0.494505494500000 }
                });

                Kcc[3] = Matrix.CreateFromArray(new double[,]
                {
                    {0.494505494500000, 0.178571428600000, 0, 0, 0, 0  },
                    {0.178571428600000, 0.494505494500000, 0, 0, 0, 0  },
                    {0, 0, 0.494505494500000, -0.178571428600000, 0, 0 },
                    {0, 0, -0.178571428600000, 0.494505494500000, 0, 0 },
                    {0, 0, 0, 0, 0.494505494500000, -0.178571428600000 },
                    {0, 0, 0, 0, -0.178571428600000, 0.494505494500000 }
                });
                return Kcc;
            }
        }

        private static Matrix MatrixKccStar => Matrix.CreateFromArray(new double[,]
        {
            {0.519272429341174, -0.0224949475741955, -0.0673726448231679, 0.0532992286325700, -0.115596477967984, -0.180005752865632, 0, 0                                   },
            {-0.0224949475741955, 0.571373230561769, 0.00636415859565785, -0.310650945646995, -0.110850590090827, -0.115596477967984, 0, 0                                   },
            {-0.0673726448231679, 0.00636415859565783, 1.13126609833817, 0, -0.323914195445245, 2.77555756156289e-17, -0.0673726448231678, -0.00636415859565788              },
            {0.0532992286325700, -0.310650945646995, 0, 1.04037205521382, 4.16333634234434e-17, -0.128818550126366, -0.0532992286325699, -0.310650945646995                  },
            {-0.115596477967984, -0.110850590090827, -0.323914195445245, 5.72458747072346e-17, 0.555107150696809, 1.38777878078145e-17, -0.115596477967984, 0.110850590090827},
            {-0.180005752865631, -0.115596477967984, 2.77555756156289e-17, -0.128818550126366, 1.38777878078145e-17, 0.360011505131263, 0.180005752865631, -0.115596477967984},
            {0, 0, -0.0673726448231679, -0.0532992286325699, -0.115596477967984, 0.180005752865631, 0.519272429341174, 0.0224949475741954                                    },
            {0, 0, -0.00636415859565789, -0.310650945646995, 0.110850590090827, -0.115596477967984, 0.0224949475741954, 0.571373230561769                                    }
        });

        private static Vector VectorDr => Vector.CreateFromArray(new double[]
        {
            0, 0, 0, 0, -3.375195492420367, -10.215251712309035, 0.418600986802971, -1.151753569240856
        });

        private static Vector VectorFcStar => Vector.CreateFromArray(new double[]
        {
            0, 0, 3.19374601989718, 2.11725876973317, 0.254044579955345, 6.55220940504195, -3.44779060118368, 1.33053183634499
        });

        private static Dictionary<int, Vector> VectorsFr
        {
            get
            {
                var fr = new Dictionary<int, Vector>();
                fr[0] = Vector.CreateZero(8);
                fr[1] = Vector.CreateZero(12);
                fr[2] = Vector.CreateZero(8);
                fr[3] = Vector.CreateZero(12);
                fr[3][11] = 10;
                return fr;
            }
        }

        private static Dictionary<int, Vector> VectorsFbc
        {
            get
            {
                var fbc = new Dictionary<int, Vector>();
                fbc[0] = Vector.CreateZero(4);
                fbc[1] = Vector.CreateZero(6);
                fbc[2] = Vector.CreateZero(4);
                fbc[3] = Vector.CreateZero(6);
                return fbc;
            }
        }

        [Fact]
        public static void TestCoarseProblem()
        {
            // Setup the model and solver
            Model model = MappingMatricesTests.CreateModel();
            Dictionary<int, INode[]> cornerNodes = MappingMatricesTests.DefineCornerNodes(model);
            var solver = new FetiDPSolver.Builder(cornerNodes).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);

            // Use the hardcoded intermediate matrices & vectors
            Dictionary<int, Vector> fbc = VectorsFbc;
            Dictionary<int, Vector> fr = VectorsFr;
            Dictionary<int, Matrix> Kcc = MatricesKcc;
            Dictionary<int, Matrix> Krc = MatricesKrc;
            Dictionary<int, Matrix> Krr = MatricesKrr;
            var factorizedKrr = new Dictionary<int, CholeskyFull>();
            for (int i = 0; i < 4; ++i) factorizedKrr[i] = Krr[i].FactorCholesky(false);

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolver).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparator)fi.GetValue(solver);

            // Calculate the coarse problem matrix and rhs
            FetiDPInterfaceProblemSolver interfaceSolver = new FetiDPInterfaceProblemSolver.Builder().Build();
            Vector globalFcStar = interfaceSolver.CreateCoarseProblemRhs(dofSeparator, factorizedKrr, Krc, fr, fbc);
            MethodInfo method = interfaceSolver.GetType().GetMethod("CreateGlobalKccStar",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Matrix globalKccStar = (Matrix)method.Invoke(interfaceSolver, new object[] { dofSeparator, factorizedKrr, Krc, Kcc });

            // Check against expected matrices
            var expectedKccStar = MatrixKccStar;
            var expectedFcStar = VectorFcStar;
            double tol = 1E-13;
            Assert.True(expectedKccStar.Equals(globalKccStar, tol));
            Assert.True(expectedFcStar.Equals(globalFcStar, tol));
        }

        [Fact]
        public static void TestDisconnectedDisplacements()
        {
            //TODO: Perhaps use the Br, Bc from the class that tests them instead of the solver.

            Model model = MappingMatricesTests.CreateModel();
            Dictionary<int, INode[]> cornerNodes = MappingMatricesTests.DefineCornerNodes(model);
            var solver = new FetiDPSolver.Builder(cornerNodes).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);

            Dictionary<int, Matrix> Krr = MatricesKrr;
            var factorizedKrr = new Dictionary<int, CholeskyFull>();
            for (int i = 0; i < 4; ++i) factorizedKrr[i] = Krr[i].FactorCholesky(false);
            var fr = VectorsFr;

            Vector dr = solver.CalcDisconnectedDisplacements(factorizedKrr, fr);
            var expectedDr = VectorDr;

            double tol = 1E-13;
            Assert.True(expectedDr.Equals(dr, tol));
        }

        [Fact]
        public static void TestFlexibilityMatrices()
        {
            // Setup the model and solver
            Model model = MappingMatricesTests.CreateModel();
            Dictionary<int, INode[]> cornerNodes = MappingMatricesTests.DefineCornerNodes(model);
            var solver = new FetiDPSolver.Builder(cornerNodes).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);

            // Use the hardcoded intermediate matrices
            Dictionary<int, Matrix> Krc = MatricesKrc;
            Dictionary<int, Matrix> Krr = MatricesKrr;
            var factorizedKrr = new Dictionary<int, CholeskyFull>();
            for (int i = 0; i < 4; ++i) factorizedKrr[i] = Krr[i].FactorCholesky(false);

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolver).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumerator)fi.GetValue(solver);
            fi = typeof(FetiDPSolver).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparator)fi.GetValue(solver);

            // Create the flexibility matrices by multiplying with identity matrices
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            int numCornerDofs = dofSeparator.NumGlobalCornerDofs;
            var flexibility = new FetiDPFlexibilityMatrix(factorizedKrr, Krc, lagrangeEnumerator, dofSeparator);
            Matrix FIrr = MultiplyWithIdentity(numLagranges, numLagranges, flexibility.MultiplyFIrr);
            Matrix FIrc = MultiplyWithIdentity(numLagranges, numLagranges, (x, y) => y.CopyFrom(flexibility.MultiplyFIrc(x)));

            // Check against expected matrices
            var expectedFIrr = Matrix.CreateFromArray(new double[,]
            {
                {3.57057200993606, -0.108283270560292, 0.338429752179871, -0.279338843056072, -0.573961878785917, -0.114111168546807, 0, 0  },
                {-0.108283270560292, 2.65633088920628, -0.234165486478537, 0.447212600200740, -0.173283461887574, -0.573961878785916, 0, 0  },
                {0.338429752179871, -0.234165486478537, 2.26748388676785, -2.77555756156289e-17, 0, 0, -0.338429752179871, -0.23416548647853},
                {-0.279338843056072, 0.447212600200740, -2.77555756156289e-17, 3.03419905760385, 0, 0, -0.279338843056072, -0.44721260020073},
                {-0.573961878785917, -0.173283461887574, 0, 0, 2.71882869786337, -2.63677968348475e-16, 0.573961878785916, -0.17328346188757},
                {-0.114111168546807, -0.573961878785916, 0, 0, -2.63677968348475e-16, 4.04692914347278, -0.114111168546807, 0.57396187878591},
                {0, 0, -0.338429752179871, -0.279338843056072, 0.573961878785916, -0.114111168546807, 3.57057200993606, 0.108283270560292   },
                {0, 0, -0.234165486478537, -0.447212600200739, -0.173283461887574, 0.573961878785916, 0.108283270560292, 2.65633088920628   }
            });

            var expectedFIrc = Matrix.CreateFromArray(new double[,]
            {
                {0.244415273977447, 0.232902352320994, 0.188150279438879, -0.367471730456911, 0.325403750731022, 0.134569378056537, 0, 0                                           },
                {-0.127173102613820, 0.0205141879116909, 0.0345524284084688, 0.0581554138518863, 0.0926206740937292, 0.0806737768451912, 0, 0                                      },
                {-0.00592361806200106, 0.0896358681318229, -5.55111512312578e-17, 0.152076488272397, 0, 0, 0.00592361806200106, 0.0896358681318228                                 },
                {0.0980092297384746, -0.270103007383575, -0.136396263304725, 0, 0, 0, 0.0980092297384746, 0.270103007383575                                                        },
                {-0.0806737768451914, -0.0926206740937293, -2.22044604925031e-16, 0.0238937950679602, -1.66533453693773e-16, 0.161347553342741, 0.0806737768451914, -0.092620674093},
                {-0.134569378056537, -0.325403750731022, 0.515640037124902, 5.55111512312578e-17, -0.246501280853069, -2.22044604925031e-16, -0.134569378056537, 0.325403750731022 },
                {0, 0, 0.188150279438879, 0.367471730456911, 0.325403750731022, -0.134569378056537, 0.244415273977447, -0.232902352320994                                          },
                {0, 0, -0.0345524284084688, 0.0581554138518863, -0.0926206740937292, 0.0806737768451912, 0.127173102613820, 0.0205141879116910                                     }
            });

            double tol = 1E-11;
            Assert.True(expectedFIrr.Equals(FIrr, tol));
            Assert.True(expectedFIrc.Equals(FIrc, tol));
        }

        [Fact]
        public static void TestInterfaceProblemCreation()
        {
            // Setup the model and solver
            Model model = MappingMatricesTests.CreateModel();
            Dictionary<int, INode[]> cornerNodes = MappingMatricesTests.DefineCornerNodes(model);
            var solver = new FetiDPSolver.Builder(cornerNodes).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);

            // Use the hardcoded intermediate matrices & vectors
            Dictionary<int, Matrix> Krc = MatricesKrc;
            Dictionary<int, Matrix> Krr = MatricesKrr;
            var factorizedKrr = new Dictionary<int, CholeskyFull>();
            for (int i = 0; i < 4; ++i) factorizedKrr[i] = Krr[i].FactorCholesky(false);
            Vector dr = VectorDr;

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolver).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumerator)fi.GetValue(solver);
            fi = typeof(FetiDPSolver).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparator)fi.GetValue(solver);

            // Hardcoded coarse problem matrix and rhs
            FetiDPInterfaceProblemSolver interfaceSolver = new FetiDPInterfaceProblemSolver.Builder().Build();
            Vector globalFcStar = VectorFcStar;
            CholeskyFull factorKccStar = MatrixKccStar.FactorCholesky(false); // It must be set as a private field using reflection.
            fi = typeof(FetiDPInterfaceProblemSolver).GetField("factorizedGlobalKccStar",
                BindingFlags.NonPublic | BindingFlags.Instance);
            fi.SetValue(interfaceSolver, factorKccStar);

            // Create the rhs vector of the interface problem 
            var flexibility = new FetiDPFlexibilityMatrix(factorizedKrr, Krc, lagrangeEnumerator, dofSeparator);
            Vector fcStar = VectorFcStar;
            MethodInfo method = interfaceSolver.GetType().GetMethod("CreateInterfaceProblemRhs",
                BindingFlags.NonPublic | BindingFlags.Instance); // reflection for the private method
            Vector interfaceRhs = (Vector)method.Invoke(interfaceSolver, new object[] { flexibility, fcStar, dr });

            // Create the matrix of the interface problem by multiplying with identity matrix
            int numLagranges = lagrangeEnumerator.NumLagrangeMultipliers;
            var interfaceMatrixImplicit = new FetiDPInterfaceProblemSolver.InterfaceProblemMatrix(flexibility, factorKccStar);
            Matrix interfaceMatrix = MultiplyWithIdentity(numLagranges, numLagranges, interfaceMatrixImplicit.Multiply); // Action<T> is contravariant!!!

            // Check against expected linear system
            double tol = 1E-13;
            var expectedInterfaceRhs = Vector.CreateFromArray(new double[]
            {
                -14.9810838729735, -5.69975426333296, -10.5434726428584, 0.244121938779135, -5.89291361392317, -13.1189445403298, 16.4060122931895, -5.93260749341458
            });
            var expectedInterfaceMatrix = Matrix.CreateFromArray(new double[,]
            {
                {4.97303228211014, 0.0303495596681853, 0.478870888249040, -0.510074267142073, -0.554638064586832, -0.639203724429228, 0.342142442028879, -0.0259960944946401        },
                {0.0303495596681853, 2.75206109564910, -0.132450092323064, 0.393630305623386, -0.137248503790822, -0.602781684493727, 0.0259960944946400, 0.0326814286491487        },
                {0.478870888249040, -0.132450092323064, 2.50773589560585, -6.12608644034637e-17, -0.00517803522230879, -1.40567757982761e-16, -0.478870888249040, -0.132450092323064},
                {-0.510074267142073, 0.393630305623386, -1.10893466892511e-16, 3.35755264805804, 1.18406464657668e-16, 0.293147921294313, -0.510074267142073, -0.393630305623386    },
                {-0.554638064586832, -0.137248503790822, -0.00517803522230876, 1.39289728666574e-16, 2.79928131129566, -1.91565648067391e-16, 0.554638064586830, -0.137248503790822 },
                {-0.639203724429228, -0.602781684493727, -2.44519425851536e-16, 0.293147921294313, -2.04349167449641e-16, 4.95200731231851, -0.639203724429227, 0.602781684493727   },
                {0.342142442028879, 0.0259960944946401, -0.478870888249040, -0.510074267142073, 0.554638064586830, -0.639203724429227, 4.97303228211014, -0.0303495596681853        },
                {-0.0259960944946401, 0.0326814286491487, -0.132450092323064, -0.393630305623386, -0.137248503790822, 0.602781684493727, -0.0303495596681853, 2.75206109564910      }
            });
            Assert.True(expectedInterfaceRhs.Equals(interfaceRhs, tol));
            Assert.True(expectedInterfaceMatrix.Equals(interfaceMatrix, tol));
        }

        [Fact]
        public static void TestInterfaceProblemSolution()
        {
            // Setup the model and solver
            Model model = MappingMatricesTests.CreateModel();
            Dictionary<int, INode[]> cornerNodes = MappingMatricesTests.DefineCornerNodes(model);
            var solver = new FetiDPSolver.Builder(cornerNodes).BuildSolver(model);
            model.ConnectDataStructures();
            solver.OrderDofs(false);

            // Use the hardcoded intermediate matrices & vectors
            Dictionary<int, Matrix> Krc = MatricesKrc;
            Dictionary<int, Matrix> Krr = MatricesKrr;
            var factorizedKrr = new Dictionary<int, CholeskyFull>();
            for (int i = 0; i < 4; ++i) factorizedKrr[i] = Krr[i].FactorCholesky(false);
            Vector dr = VectorDr;

            // Access private fields of FetiDPSolver
            FieldInfo fi = typeof(FetiDPSolver).GetField("lagrangeEnumerator", BindingFlags.NonPublic | BindingFlags.Instance);
            var lagrangeEnumerator = (FetiDPLagrangeMultipliersEnumerator)fi.GetValue(solver);
            fi = typeof(FetiDPSolver).GetField("dofSeparator", BindingFlags.NonPublic | BindingFlags.Instance);
            var dofSeparator = (FetiDPDofSeparator)fi.GetValue(solver);

            // Hardcoded coarse problem matrix and rhs
            FetiDPInterfaceProblemSolver interfaceSolver = new FetiDPInterfaceProblemSolver.Builder().Build();
            Vector globalFcStar = VectorFcStar;
            CholeskyFull factorKccStar = MatrixKccStar.FactorCholesky(false); // It must be set as a private field using reflection.
            fi = typeof(FetiDPInterfaceProblemSolver).GetField("factorizedGlobalKccStar", 
                BindingFlags.NonPublic | BindingFlags.Instance);
            fi.SetValue(interfaceSolver, factorKccStar);

            // Dirichlet preconditioner
            var precondFactory = new DirichletPreconditioner.Factory();
            var repackagedKrr = new Dictionary<int, IMatrixView>();
            foreach (var idMatrixPair in Krr) repackagedKrr[idMatrixPair.Key] = idMatrixPair.Value;
            var stiffnessDistribution = new HomogeneousStiffnessDistribution(model, dofSeparator);
            IFetiPreconditioner preconditioner = 
                precondFactory.CreatePreconditioner(stiffnessDistribution, dofSeparator, lagrangeEnumerator, repackagedKrr);

            // Solve the interface problem
            var flexibility = new FetiDPFlexibilityMatrix(factorizedKrr, Krc, lagrangeEnumerator, dofSeparator);
            double globalForcesNorm = 10.0; // sqrt(0^2 + 0^2 + ... + 0^2 + 10^2)
            var logger = new DualSolverLogger();
            (Vector lagranges, Vector uc) = 
                interfaceSolver.SolveInterfaceProblem(flexibility, preconditioner, globalFcStar, dr, globalForcesNorm, logger);

            // Check against expected solution
            double tol = 1E-7;
            Vector expectedLagranges = Vector.CreateFromArray(new double[]
            {
                -3.67505611805653, -3.06916047739931, -3.12635180105707, 0.427127980701075, -3.73923329344533, -2.87580179407164, 3.34727977833535, -1.76301688321532
            });
            Vector expectedUc = Vector.CreateFromArray(new double[]
            {
                21.1181096194325, 27.2052778266603, 1.63160365361360, 25.0125476374046, 1.69318450304898, 61.8894615542161, -24.7267309921556, 26.3640977652349
            });
            Assert.True(expectedLagranges.Equals(lagranges, tol));
            Assert.True(expectedUc.Equals(uc, tol));
        }

        [Fact]
        public static void TestSolver()
        {
            // Setup the model and solver
            Model model = MappingMatricesTests.CreateModel();
            Dictionary<int, INode[]> cornerNodes = MappingMatricesTests.DefineCornerNodes(model);
            var solver = new FetiDPSolver.Builder(cornerNodes).BuildSolver(model);
            var problem = new ProblemStructural(model, solver);
            var linearAnalyzer = new LinearAnalyzer(model, solver, problem);
            var staticAnalyzer = new StaticAnalyzer(model, solver, problem, linearAnalyzer);

            // Run the analysis
            staticAnalyzer.Initialize();
            staticAnalyzer.Solve();

            // Gather the global displacements
            var sudomainDisplacements = new Dictionary<int, IVectorView>();
            foreach (var ls in solver.LinearSystems) sudomainDisplacements[ls.Key] = ls.Value.Solution;
            Vector globalU = solver.GatherGlobalDisplacements(sudomainDisplacements);

            // Check against expected solution
            double tol = 1E-7;
            var globalUExpected = Vector.CreateFromArray(new double[]
            {
                13.3258563908201, 12.3999624809163, 21.1181095863809, 27.2052777811441, 24.3525812415758,
                43.2777053704649, 24.8992347210378, 57.3521080292628, 4.74521148903743, 9.87352108397423,
                9.37569840211612, 25.9120840082139, 11.8910608093164, 43.4314599456699, 12.2652060584230,
                57.9466725072280, 0.450346260334126, 9.02020634682474, 1.63160365355026, 25.0125475922504,
                2.58948267402381, 45.0651412625480, 1.69318450300533, 61.8894614604312, -4.68849826343688,
                8.90417219731433, -8.76400355420594, 24.5661224138922, -9.40948533633272, 47.1084814579881,
                -11.2141368968962, 73.2559168929990, -14.0271645568764, 11.3572597884005, -24.7267309592324,
                26.3640977197317, -34.3702668180117, 46.6307017985724, -42.8927307907656, 96.0971764416081
            });
            Assert.True(globalUExpected.Equals(globalU, tol));
        }

        private static Matrix MultiplyWithIdentity(int numRows, int numCols, Action<Vector, Vector> matrixVectorMultiplication)
        {
            var result = Matrix.CreateZero(numRows, numCols);
            for (int j = 0; j < numCols; ++j)
            {
                var lhs = Vector.CreateZero(numCols);
                lhs[j] = 1.0;
                var rhs = Vector.CreateZero(numRows);
                matrixVectorMultiplication(lhs, rhs);
                result.SetSubcolumn(j, rhs);
            }
            return result;
        }
    }
}
