using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses
{
    /// <summary>
    /// Mesh and model creation supportive class, that can be used for 3d or 2d rve problems or structured meshing of parts in general
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class FEMMeshBuilder   
    {
        public static int[,] topologia_shell_coh(int elements, int elem1, int elem2, object komvoi_8)
        {
            int elem;
            int[,] t_shell = new int[elements, 8];
            for (int nrow = 0; nrow < elem1; nrow++)
            {
                for (int nline = 0; nline < elem2; nline++)
                {
                    elem = (nrow + 1 - 1) * elem2 + nline + 1;//nrow+ 1 nline+1 einai zero based 
                    t_shell[elem - 1, -1 + 1] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 8] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 4] = (nrow + 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;

                    t_shell[elem - 1, -1 + 5] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 2;
                    t_shell[elem - 1, -1 + 7] = (nrow + 1 - 1) * (3 * elem2 + 2) + 2 * elem2 + 1 + (nline + 1 - 1) * 1 + 1;

                    t_shell[elem - 1, -1 + 2] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 3;
                    t_shell[elem - 1, -1 + 6] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 2;
                    t_shell[elem - 1, -1 + 3] = (nrow + 1 - 1) * (3 * elem2 + 2) + (nline + 1 - 1) * 2 + 1;
                }
            }
            return t_shell;

        }

        public static int Topol_rve(int h1, int h2, int h3, int hexa1, int hexa2, int hexa3, int kuvos, int endiam_plaka, int katw_plaka)
        {
            int arith;
            if (h3 == 1)
            { arith = h1 + (h2 - 1) * (hexa1 + 1) + kuvos; }
            else
            {
                if (h3 == hexa3 + 1)
                { arith = hexa3 * (hexa1 + 1) * (hexa2 + 1) + h1 + (h2 - 1) * (hexa1 + 1); }
                else
                {
                    if (h2 == 1)
                    { arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + h1; }
                    else
                    {
                        if (h2 == hexa2 + 1)
                        { arith = (h3 - 2) * endiam_plaka + kuvos + katw_plaka + (hexa1 + 1) + 2 * (hexa2 - 1) + h1; }
                        else
                        {
                            if (h1 == 1)
                            { arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 1; }
                            else
                            {
                                if (h1 == hexa1 + 1)
                                { arith = kuvos + katw_plaka + (h3 - 2) * endiam_plaka + (hexa1 + 1) + (h2 - 2) * 2 + 2; }
                                else
                                { arith = (h1 - 1) + (h2 - 2) * (hexa1 - 1) + (h3 - 2) * (hexa1 - 1) * (hexa2 - 1); }
                            }
                        }
                    }

                }
            }
            return arith;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParameters(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 95, //150, // diastaseis
                L02 = 95, //150,
                L03 = 95, //40,
                hexa1 = discr1 * subdiscr1,// diakritopoihsh
                hexa2 = discr1 * subdiscr1,
                hexa3 = discr3
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = discr1_shell * subdiscr1_shell,
                elem2 = discr1_shell * subdiscr1_shell,
                L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
                L2 = 50,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.05,// Gpa
                D_o_1 = 0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceKanonikhGewmetriaRveExampleParametersStiffCase(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 95, //150, // diastaseis
                L02 = 95, //150,
                L03 = 95, //40,
                hexa1 = discr1 * subdiscr1,// diakritopoihsh
                hexa2 = discr1 * subdiscr1,
                hexa3 = discr3
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = discr1_shell * subdiscr1_shell,
                elem2 = discr1_shell * subdiscr1_shell,
                L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
                L2 = 50,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.20, //0.05,  // 1Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.25, //0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.20, //0.05,// Gpa
                D_o_1 = 0.25, //0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static Tuple<rveMatrixParameters, grapheneSheetParameters> GetReferenceRveExampleParameters(int subdiscr1, int discr1, int discr3, int subdiscr1_shell, int discr1_shell)
        {
            rveMatrixParameters mp;
            mp = new rveMatrixParameters()
            {
                E_disp = 3.5, //Gpa
                ni_disp = 0.4, // stather Poisson
                L01 = 100, //150, // diastaseis
                L02 = 100, //150,
                L03 = 100, //40,
                hexa1 = discr1 * subdiscr1,// diakritopoihsh
                hexa2 = discr1 * subdiscr1,
                hexa3 = discr3
            };

            grapheneSheetParameters gp;
            gp = new grapheneSheetParameters()
            {
                // parametroi shell
                E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
                ni_shell = 0.0607, // stathera poisson
                elem1 = discr1_shell * subdiscr1_shell,
                elem2 = discr1_shell * subdiscr1_shell,
                L1 = 50,// nm  // DIORTHOSI 2 graphene sheets
                L2 = 50,// nm
                L3 = 112.5096153846, // nm
                a1_shell = 0, // nm
                tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

                //parametroi cohesive epifaneias
                T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
                D_o_3 = 0.5, // nm
                D_f_3 = 4, // nm
                T_o_1 = 0.05,// Gpa
                D_o_1 = 0.5, // nm
                D_f_1 = 4, // nm
                n_curve = 1.4
            };

            Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
            return gpmp;
        }

        public static void HexaElementsOnlyRVEwithRenumbering_forMS(Model model, rveMatrixParameters mp, double[,] Dq, string renumberingVectorPath, Dictionary<int, Node> boundaryNodes)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Rve Matrix
            double E_disp = mp.E_disp; //Gpa
            double ni_disp = mp.ni_disp; // stather Poisson
            double L01 = mp.L01; // diastaseis
            double L02 = mp.L02;
            double L03 = mp.L03;
            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            // Perioxh parametroi Rve Matrix ews edw


            // Perioxh Gewmetria shmeiwn
            int nodeCounter = 0;

            int nodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int h1 = 0; h1 < hexa1 + 1; h1++)
            {
                for (int h2 = 0; h2 < hexa2 + 1; h2++)
                {
                    for (int h3 = 0; h3 < hexa3 + 1; h3++)
                    {
                        nodeID = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                        nodeCoordY =-0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                        nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                        model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y: nodeCoordY, z: nodeCoordZ ));
                        nodeCounter++;
                    }
                }
            }
            // Perioxh Gewmetria shmeiwn ews edw

            //Perioxh Eisagwgh elements
            int elementCounter = 0;
            int subdomainID = 1;

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };
            Element e1;
            int ElementID;
            int[] globalNodeIDforlocalNode_i = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
                        globalNodeIDforlocalNode_i[0] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[1] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[2] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[3] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[4] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[5] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[6] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[7] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

                        e1 = new Element()
                        {
                            ID = ElementID,
                            ElementType = new Hexa8NonLinear(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                        };

                        for (int j = 0; j < 8; j++)
                        {
                            e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
                        }
                        model.ElementsDictionary.Add(e1.ID, e1);
                        model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
                        elementCounter++;
                    }
                }
            }
            //Perioxh Eisagwgh elements

            //Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
            // change one tuple value
            //nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
            // get one tuple value
            //elementCounter = nodeElementCounters.Item2;            
            //return nodeElementCounters;

            int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int f_komvoi_rve = kuvos;
            int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            int komvos;
            Dq = new double[9, 3 * p_komvoi_rve];
            for (int j = 0; j < p_komvoi_rve; j++)
            {
                komvos = renumbering.GetNewNodeNumbering(f_komvoi_rve + j + 1);
                Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
                Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
                Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
                Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
                Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
                Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
                Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
                Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
                Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
                boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
            }
        }

        public static void AddGrapheneSheet_with_o_x_Input_withRenumbering(Model model, grapheneSheetParameters gp, double[] ekk_xyz, o_x_parameters o_x_parameters, string renumberingVectorPath, string o_xsunol_input_path)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            double[] o_xsunol;
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            //o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz, o_x_parameters);
            o_xsunol = PrintUtilities.ReadVector(o_xsunol_input_path);

            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY = o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ));
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            var material2 = new ShellElasticMaterial3D()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK = 5 / 6,
            };

            int elements = elem1 * elem2;
            int fdof_8 = 5 * (elem1 * (3 * elem2 + 2) + 2 * elem2 + 1);
            int komvoi_8 = fdof_8 / 5;
            int[,] t_shell;
            t_shell = FEMMeshBuilder.topologia_shell_coh(elements, elem1, elem2, komvoi_8); // ta stoixeia tou einai 1 based to idio einai 0 based

            double[] Tk_vec = new double[8];
            double[][] VH = new double[8][];
            int[] midsurfaceNodeIDforlocalShellNode_i = new int[8];
            Element e2;
            int ElementID;
            int subdomainID = 1;

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }

            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalShellNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    //
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) //ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
                    {
                        //oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY =o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ));
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BenzeggaghKenaneCohesiveMaterial material3 = new Materials.BenzeggaghKenaneCohesiveMaterial()
            {
                T_o_3 = T_o_3,
                D_o_3 = D_o_3,
                D_f_3 = D_f_3,
                T_o_1 = T_o_1,
                D_o_1 = D_o_1,
                D_f_1 = D_f_1,
                n_curve = n_curve,
            };

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].Elements.Add( e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY =o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ ));
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].Elements.Add( e2);
                eswterikosElementCounter++;
            }
            // orismos elements anw strwshs ews edw

        }

        public static Dictionary<Node, IList<IDofType>> GetConstraintsOfDegenerateRVEForNonSingularStiffnessMatrix_withRenumbering(Model model, int hexa1, int hexa2, int hexa3, string renumberingVectorPath)
        {
            //Origin : RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix_withRenumbering()
            //modifications: return type and nodes and dofs to be constrained

            Dictionary<Node, IList<IDofType>> RigidBodyNodeConstraints = new Dictionary<Node, IList<IDofType>>();

            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);
            int nodeID;

            nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
            RigidBodyNodeConstraints.Add(model.NodesDictionary[nodeID], new List<IDofType>() { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });


            nodeID = renumbering.GetNewNodeNumbering(Topol_rve(hexa1 + 1, 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
            RigidBodyNodeConstraints.Add(model.NodesDictionary[nodeID], new List<IDofType>() { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });

            nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, hexa2 + 1, 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
            RigidBodyNodeConstraints.Add(model.NodesDictionary[nodeID], new List<IDofType>() { StructuralDof.TranslationX, StructuralDof.TranslationY, StructuralDof.TranslationZ });

            //nodeID = renumbering.GetNewNodeNumbering(Topol_rve(1, 1, hexa3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

            return RigidBodyNodeConstraints;

        }

        public static void AddGrapheneSheet_with_o_x_Input_withRenumberingBondSlip(Model model, grapheneSheetParameters gp, double[] ekk_xyz, o_x_parameters o_x_parameters, string renumberingVectorPath, string o_xsunol_input_path)
        {
            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Graphene sheet
            // parametroi shell
            double E_shell = gp.E_shell; // GPa = 1000Mpa = 1000N / mm2
            double ni_shell = gp.ni_shell; // stathera poisson
            int elem1 = gp.elem1;
            int elem2 = gp.elem2;
            double L1 = gp.L1;// nm
            double L2 = gp.L2;// nm
            double L3 = gp.L3; // nm
            double a1_shell = gp.a1_shell; // nm
            double tk = gp.tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            double T_o_3 = gp.T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            double D_o_3 = gp.D_o_3; // nm
            double D_f_3 = gp.D_f_3; // nm

            double T_o_1 = gp.T_o_1;// Gpa
            double D_o_1 = gp.D_o_1; // nm
            double D_f_1 = gp.D_f_1; // nm

            double n_curve = gp.n_curve;
            // Perioxh parametroi Graphene sheet ews edw


            int eswterikosNodeCounter = 0;
            int eswterikosElementCounter = 0;
            int PreviousElementsNumberValue = model.ElementsDictionary.Count();
            int PreviousNodesNumberValue = model.NodesDictionary.Count();


            // Perioxh gewmetrias (orismos nodes) meshs epifaneias
            int new_rows = 2 * elem1 + 1;
            int new_lines = 2 * elem2 + 1;
            double[] o_xsunol;
            int NodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;

            //o_xsunol = ox_sunol_Builder_ekk_with_o_x_parameters(new_rows, new_lines, L1, L2, elem1, elem2, a1_shell, ekk_xyz, o_x_parameters);
            o_xsunol = PrintUtilities.ReadVector(o_xsunol_input_path);

            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0];
                nodeCoordY =o_xsunol[6 * nNode + 1];
                nodeCoordZ = o_xsunol[6 * nNode + 2];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ ));
                eswterikosNodeCounter++;
            }
            int arithmosShmeiwnShellMidsurface = eswterikosNodeCounter;
            // perioxh gewmetrias meshs epifaneias ews edw


            // perioxh orismou shell elements
            var material2 = new ShellElasticMaterial3D()
            {
                YoungModulus = E_shell,
                PoissonRatio = ni_shell,
                ShearCorrectionCoefficientK = 5 / 6,
            };

            int elements = elem1 * elem2;
            int fdof_8 = 5 * (elem1 * (3 * elem2 + 2) + 2 * elem2 + 1);
            int komvoi_8 = fdof_8 / 5;
            int[,] t_shell;
            t_shell = FEMMeshBuilder.topologia_shell_coh(elements, elem1, elem2, komvoi_8); // ta stoixeia tou einai 1 based to idio einai 0 based

            double[] Tk_vec = new double[8];
            double[][] VH = new double[8][];
            int[] midsurfaceNodeIDforlocalShellNode_i = new int[8];
            Element e2;
            int ElementID;
            int subdomainID = 1;

            for (int j = 0; j < 8; j++) // paxos idio gia ola telements
            {
                Tk_vec[j] = tk;
            }

            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalShellNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    //
                    ElementType = new Shell8NonLinear(material2, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) //ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
                    {
                        //oVn_i= new double[][] { new double [] {ElementID, ElementID }, new double [] { ElementID, ElementID } },
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalShellNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalShellNode_i[j1] + PreviousNodesNumberValue)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].Elements.Add( e2);
                eswterikosElementCounter++;
            }
            int arithmosShellElements = eswterikosElementCounter;
            // perioxh orismou shell elements ews edw

            // orismos shmeiwn katw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] - 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] - 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] - 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ ));
                eswterikosNodeCounter++;
            }
            //

            //orismos elements katw strwshs
            BondSlipCohMat material3 = new Materials.BondSlipCohMat(T_o_1,D_o_1,0.1,T_o_3,D_o_3,new double[2],new double[2],1e-10);            

            int[] midsurfaceNodeIDforlocalCohesiveNode_i = new int[8];
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 0,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].Elements.Add( e2);
                eswterikosElementCounter++;
            }
            // orismos elements katw strwshs ews edw

            // orismos shmeiwn anw strwshs
            for (int nNode = 0; nNode < o_xsunol.GetLength(0) / 6; nNode++) //nNode einai zero based
            {
                NodeID = renumbering.GetNewNodeNumbering(eswterikosNodeCounter + PreviousNodesNumberValue + 1);
                nodeCoordX = o_xsunol[6 * nNode + 0] + 0.5 * tk * o_xsunol[6 * nNode + 3];
                nodeCoordY = o_xsunol[6 * nNode + 1] + 0.5 * tk * o_xsunol[6 * nNode + 4];
                nodeCoordZ = o_xsunol[6 * nNode + 2] + 0.5 * tk * o_xsunol[6 * nNode + 5];

                model.NodesDictionary.Add(NodeID, new Node(id: NodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ ));
                eswterikosNodeCounter++;
            }
            //
            //orismos elements anw strwshs 
            //TODO material3 = new Materials.BondSlipCohMat(T_o_1, D_o_1, 0.1, T_o_3, D_o_3, new double[2], new double[2], 1e-10);
            for (int nElement = 0; nElement < elements; nElement++)
            {
                ElementID = eswterikosElementCounter + PreviousElementsNumberValue + 1;
                // ta dianusmata katefthunshs allazoun analoga to element 
                for (int j1 = 0; j1 < 8; j1++)
                {
                    midsurfaceNodeIDforlocalCohesiveNode_i[j1] = t_shell[nElement, j1]; // periexei NOT zero based 
                    VH[j1] = new double[3];
                    VH[j1][0] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 3];
                    VH[j1][1] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 4];
                    VH[j1][2] = o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[j1] - 1) + 5];
                }

                e2 = new Element()
                {
                    ID = ElementID,
                    ElementType = new CohesiveShell8ToHexa20(material3, GaussLegendre2D.GetQuadratureWithOrder(3, 3)) //ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_1(material3, 3, 3)
                    {
                        oVn_i = new double[][] { new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[0] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[1] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[2] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[3] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[4] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[5] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[6] - 1) + 5] },
                                                 new double[] { o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 3], o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 4],o_xsunol[6 * (midsurfaceNodeIDforlocalCohesiveNode_i[7] - 1) + 5] },},
                        tk = Tk_vec,
                        ShellElementSide = 1,
                    }
                };
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue), model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue)]);
                }
                for (int j1 = 0; j1 < 8; j1++)
                {
                    e2.NodesDictionary.Add(renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface),
                        model.NodesDictionary[renumbering.GetNewNodeNumbering(midsurfaceNodeIDforlocalCohesiveNode_i[j1] + PreviousNodesNumberValue + 2 * arithmosShmeiwnShellMidsurface)]);
                }
                model.ElementsDictionary.Add(e2.ID, e2);
                model.SubdomainsDictionary[subdomainID].Elements.Add( e2);
                eswterikosElementCounter++;
            }

            // orismos elements anw strwshs ews edw

        }

        public static void LinearHexaElementsOnlyRVEwithRenumbering_forMS(Model model, rveMatrixParameters mp, double[,] Dq, string renumberingVectorPath, Dictionary<int, Node> boundaryNodes)
        {
            //COPY apo FEMMeshBuilder.HexaElementsOnlyRVEwithRenumbering_forMS()
            //modifications grammika elements kai artihmisi nodes pou afta xreiazontai

            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Rve Matrix
            double E_disp = mp.E_disp; //Gpa
            double ni_disp = mp.ni_disp; // stather Poisson
            double L01 = mp.L01; // diastaseis
            double L02 = mp.L02;
            double L03 = mp.L03;
            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            // Perioxh parametroi Rve Matrix ews edw


            // Perioxh Gewmetria shmeiwn
            int nodeCounter = 0;

            int nodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int h1 = 0; h1 < hexa1 + 1; h1++)
            {
                for (int h2 = 0; h2 < hexa2 + 1; h2++)
                {
                    for (int h3 = 0; h3 < hexa3 + 1; h3++)
                    {
                        nodeID = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                        nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                        nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                        model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ ));
                        nodeCounter++;
                    }
                }
            }
            // Perioxh Gewmetria shmeiwn ews edw

            //Perioxh Eisagwgh elements
            int elementCounter = 0;
            int subdomainID = 1;

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };
            Element e1;
            int ElementID;
            int[] globalNodeIDforlocalNode_i = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
                        globalNodeIDforlocalNode_i[6] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[7] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[4] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[5] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[2] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[3] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[0] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[1] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

                        e1 = new Element()
                        {
                            ID = ElementID,
                            ElementType = new Hexa8Fixed(material1),//, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) 
                        };

                        for (int j = 0; j < 8; j++)
                        {
                            e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
                        }
                        model.ElementsDictionary.Add(e1.ID, e1);
                        model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
                        elementCounter++;
                    }
                }
            }
            //Perioxh Eisagwgh elements

            //Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
            // change one tuple value
            //nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
            // get one tuple value
            //elementCounter = nodeElementCounters.Item2;            
            //return nodeElementCounters;

            //TODO: afta den xreiazontai poia (Dq klp)

            int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            int f_komvoi_rve = kuvos;
            int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            int komvos;
            Dq = new double[9, 3 * p_komvoi_rve];
            for (int j = 0; j < p_komvoi_rve; j++)
            {
                komvos = renumbering.GetNewNodeNumbering(f_komvoi_rve + j + 1);
                Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
                Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
                Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
                Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
                Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
                Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
                Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
                Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
                Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
                boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
            }
        }

        public static void LinearHexaElementsOnlyRVEwithRenumbering_forMS_PeripheralNodes(Model model, rveMatrixParameters mp, double[,] Dq, string renumberingVectorPath, Dictionary<int, Node> boundaryNodes)
        {
            //COPY apo FEMMeshBuilder.LinearHexaElementsOnlyRVEwithRenumbering_forMS()
            //modifications boundary nodes mono ta peripheral 

            // Perioxh renumbering initialization 
            renumbering renumbering = new renumbering(PrintUtilities.ReadIntVector(renumberingVectorPath));
            // perioxh renumbering initialization ews edw 

            // Perioxh parametroi Rve Matrix
            double E_disp = mp.E_disp; //Gpa
            double ni_disp = mp.ni_disp; // stather Poisson
            double L01 = mp.L01; // diastaseis
            double L02 = mp.L02;
            double L03 = mp.L03;
            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            // Perioxh parametroi Rve Matrix ews edw


            // Perioxh Gewmetria shmeiwn
            int nodeCounter = 0;

            int nodeID;
            double nodeCoordX;
            double nodeCoordY;
            double nodeCoordZ;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            for (int h1 = 0; h1 < hexa1 + 1; h1++)
            {
                for (int h2 = 0; h2 < hexa2 + 1; h2++)
                {
                    for (int h3 = 0; h3 < hexa3 + 1; h3++)
                    {
                        nodeID = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka)); // h1+1 dioti h1 einai zero based
                        nodeCoordX = -0.5 * L01 + (h1 + 1 - 1) * (L01 / hexa1);  // h1+1 dioti h1 einai zero based
                        nodeCoordY = -0.5 * L02 + (h2 + 1 - 1) * (L02 / hexa2);
                        nodeCoordZ = -0.5 * L03 + (h3 + 1 - 1) * (L03 / hexa3);

                        model.NodesDictionary.Add(nodeID, new Node(id: nodeID, x: nodeCoordX, y:  nodeCoordY, z: nodeCoordZ ));
                        nodeCounter++;
                    }
                }
            }
            // Perioxh Gewmetria shmeiwn ews edw

            //Perioxh Eisagwgh elements
            int elementCounter = 0;
            int subdomainID = 1;

            ElasticMaterial3D material1 = new ElasticMaterial3D()
            {
                YoungModulus = E_disp,
                PoissonRatio = ni_disp,
            };
            Element e1;
            int ElementID;
            int[] globalNodeIDforlocalNode_i = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based
                        globalNodeIDforlocalNode_i[6] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[7] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[4] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[5] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[2] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[3] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[0] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                        globalNodeIDforlocalNode_i[1] = renumbering.GetNewNodeNumbering(Topol_rve(h1 + 1 + 1, h2 + 1, h3 + 1, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));

                        e1 = new Element()
                        {
                            ID = ElementID,
                            ElementType = new Hexa8Fixed(material1),//, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) 
                        };

                        for (int j = 0; j < 8; j++)
                        {
                            e1.NodesDictionary.Add(globalNodeIDforlocalNode_i[j], model.NodesDictionary[globalNodeIDforlocalNode_i[j]]);
                        }
                        model.ElementsDictionary.Add(e1.ID, e1);
                        model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
                        elementCounter++;
                    }
                }
            }
            //Perioxh Eisagwgh elements

            //Tuple<int, int> nodeElementCounters = new Tuple<int, int>(nodeCounter, elementCounter);
            // change one tuple value
            //nodeElementCounters = new Tuple<int, int>(nodeElementCounters.Item1, elementCounter);
            // get one tuple value
            //elementCounter = nodeElementCounters.Item2;            
            //return nodeElementCounters;

            //TODO: afta den xreiazontai poia (Dq klp)

            //int komvoi_rve = (hexa1 + 1) * (hexa2 + 1) * (hexa3 + 1);
            //int f_komvoi_rve = kuvos;
            //int p_komvoi_rve = komvoi_rve - f_komvoi_rve;
            //int komvos;
            //Dq = new double[9, 3 * p_komvoi_rve];
            //for (int j = 0; j < p_komvoi_rve; j++)
            //{
            //    komvos = renumbering.GetNewNodeNumbering(f_komvoi_rve + j + 1);
            //    Dq[0, 3 * j + 0] = model.NodesDictionary[komvos].X;
            //    Dq[1, 3 * j + 1] = model.NodesDictionary[komvos].Y;
            //    Dq[2, 3 * j + 2] = model.NodesDictionary[komvos].Z;
            //    Dq[3, 3 * j + 0] = model.NodesDictionary[komvos].Y;
            //    Dq[4, 3 * j + 1] = model.NodesDictionary[komvos].Z;
            //    Dq[5, 3 * j + 2] = model.NodesDictionary[komvos].X;
            //    Dq[6, 3 * j + 0] = model.NodesDictionary[komvos].Z;
            //    Dq[7, 3 * j + 1] = model.NodesDictionary[komvos].X;
            //    Dq[8, 3 * j + 2] = model.NodesDictionary[komvos].Y;
            //    boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
            //}

            int i_2 = 1;
            for (int i1=1; i1<hexa1+2;i1++)
            {
                for (int i3 = 1; i3 < hexa3 + 2; i3++)
                {
                    int komvos = renumbering.GetNewNodeNumbering(Topol_rve(i1, i_2, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                    boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
                }
            }

            i_2 = hexa2+1;
            for (int i1 = 1; i1 < hexa1 + 2; i1++)
            {
                for (int i3 = 1; i3 < hexa3 + 2; i3++)
                {
                    int komvos = renumbering.GetNewNodeNumbering(Topol_rve(i1, i_2, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                    boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
                }
            }

            int i_1 =  1;
            for (int i2 = 2; i2 < hexa2 + 1; i2++)
            {
                for (int i3 = 1; i3 < hexa3 + 2; i3++)
                {
                    int komvos = renumbering.GetNewNodeNumbering(Topol_rve(i_1, i2, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                    boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
                }
            }

            i_1 = hexa1+1;
            for (int i2 = 2; i2 < hexa2 + 1; i2++)
            {
                for (int i3 = 1; i3 < hexa3 + 2; i3++)
                {
                    int komvos = renumbering.GetNewNodeNumbering(Topol_rve(i_1, i2, i3, hexa1, hexa2, hexa3, kuvos, endiam_plaka, katw_plaka));
                    boundaryNodes.Add(komvos, model.NodesDictionary[komvos]);
                }
            }
        }

        public static IEnumerable<Element> GetHostGroupForCohesiveElement(Element cohesive, rveMatrixParameters mp, Model model, string renumberingVectorPath)
        {            

            int hexa1 = mp.hexa1;// diakritopoihsh
            int hexa2 = mp.hexa2;
            int hexa3 = mp.hexa3;
            int kuvos = (hexa1 - 1) * (hexa2 - 1) * (hexa3 - 1);
            int endiam_plaka = 2 * (hexa1 + 1) + 2 * (hexa2 - 1);
            int katw_plaka = (hexa1 + 1) * (hexa2 + 1);

            IList<Element> possibleHosts = new List<Element>(8);

            double tol = 3e-10; //apo hexa8.GetNaturalCoordinates

            //RVE data
            double L_1 = mp.L01 / ((double)mp.hexa1);
            double L_2 = mp.L02 / ((double)mp.hexa2);
            double L_3 = mp.L03 / ((double)mp.hexa3);

            var cohNodes = cohesive.NodesDictionary.Values.ToList();
            for (int i1 = 8; i1 < 16; i1++)
            {
                Node embNode = cohNodes[i1];

                // kata to dhmiourgia_sundesmologias_embed_t.m
                double x_ek = embNode.X + 0.5 * mp.L01;
                double y_ek = embNode.Y + 0.5 * mp.L02;
                double z_ek = embNode.Z + 0.5 * mp.L03;

                var nhexa_i = new int[3][];

                var x_div = Math.Truncate(x_ek / L_1); // 'div gia double'                
                var x_pres1 = x_ek - (double)x_div * L_1;
                var x_pres2 = -x_ek + ((double)x_div + 1) * L_1;

                if (x_pres1 < tol)
                {
                    nhexa_i[0] = new int[2] { (int)x_div + 1, (int)x_div };
                }
                else
                {
                    if (x_pres2 < tol)
                    {
                        nhexa_i[0] = new int[2] { (int)x_div + 1, (int)x_div + 2 };
                    }
                    else
                    {
                        nhexa_i[0] = new int[1] { (int)x_div + 1 };
                    }
                }

                var y_div = Math.Truncate(y_ek / L_2); // 'div gia double'                
                var y_pres1 = y_ek - (double)y_div * L_2;
                var y_pres2 = -y_ek + ((double)y_div + 1) * L_2;

                if (y_pres1 < tol)
                {
                    nhexa_i[1] = new int[2] { (int)y_div + 1, (int)y_div };
                }
                else
                {
                    if (y_pres2 < tol)
                    {
                        nhexa_i[1] = new int[2] { (int)y_div + 1, (int)y_div + 2 };
                    }
                    else
                    {
                        nhexa_i[1] = new int[1] { (int)y_div + 1 };
                    }
                }


                var z_div = Math.Truncate(z_ek / L_3); // 'div gia double'                
                var z_pres1 = z_ek - (double)z_div * L_3;
                var z_pres2 = -z_ek + ((double)z_div + 1) * L_3;

                if (z_pres1 < tol)
                {
                    nhexa_i[2] = new int[2] { (int)z_div + 1, (int)z_div };
                }
                else
                {
                    if (z_pres2 < tol)
                    {
                        nhexa_i[2] = new int[2] { (int)z_div + 1, (int)z_div + 2 };
                    }
                    else
                    {
                        nhexa_i[2] = new int[1] { (int)z_div + 1 };
                    }
                }

                //Correction
                //nhexa_i[0] = new int[3] { (int)x_div , (int)x_div + 1,(int)x_div + 2 }; nhexa_i[1] = new int[3] { (int)y_div  , (int)y_div + 1,(int)y_div + 2 }; nhexa_i[2] = new int[3] { (int)z_div , (int)z_div + 1,(int)z_div+ 2 };

                for (int j1=0; j1<nhexa_i[0].GetLength(0);j1++)
                {
                    for (int j2 = 0; j2 < nhexa_i[1].GetLength(0); j2++)
                    {
                        for (int j3 = 0; j3 < nhexa_i[2].GetLength(0); j3++)
                        {
                            int possibleHostId= nhexa_i[0][j1] + (nhexa_i[1][j2] - 1) * hexa1 + (nhexa_i[2][j3] - 1) * (hexa1) * hexa2;

                            if (!possibleHosts.Contains(model.ElementsDictionary[possibleHostId]))
                            {
                                possibleHosts.Add(model.ElementsDictionary[possibleHostId]);
                            }
                        }
                    }
                }

            }

            List<int> possibleHostIds = new List<int>(possibleHosts.Count());
            for (int i1 = 0; i1 < possibleHosts.Count; i1++)
            {
                possibleHostIds.Add(possibleHosts.ElementAt(i1).ID);
            }
            possibleHostIds.Sort();
            possibleHosts = new List<Element>(possibleHosts.Count());
            for (int k1 = 0; k1 < possibleHostIds.Count; k1++)
            {
                possibleHosts.Add(model.ElementsDictionary[possibleHostIds[k1]]);
            }
            return possibleHosts;
        
        }
    }
}
