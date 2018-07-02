using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA
{
    public class ModelCreator
    {
        
        public Model Model { get; private set; }
        public int NumberOfDimensions { get; set; }
        public double Thickness { get; set; }
        public int NumberOfControlPoints { get; set; }
        public int NumberOfPatches { get; set; }
        public IFiniteElementMaterial Material {get;set;}
        public Dictionary<int, int> DegreeKsiDictionary=new Dictionary<int, int>();
        public Dictionary<int, int> DegreeHetaDictionary = new Dictionary<int, int>();
        public Dictionary<int, int> DegreeZetaDictionary = new Dictionary<int, int>();
        public Dictionary<int, double[]> KnotValueVectorsKsiDictionary = new Dictionary<int, double[]>();
        public Dictionary<int, double[]> KnotValueVectorsHetaDictionary = new Dictionary<int, double[]>();
        public Dictionary<int, double[]> KnotValueVectorsZetaDictionary = new Dictionary<int, double[]>();
        public Dictionary<int, int[]> ControlPointIDsDictionary = new Dictionary<int, int[]>();
        public Dictionary<int, int> NumberOfControlPointsKsiDictionary = new Dictionary<int, int>();
        public Dictionary<int, int> NumberOfControlPointsHetaDictionary = new Dictionary<int, int>();
        public Dictionary<int, int> NumberOfControlPointsZetaDictionary = new Dictionary<int, int>();
        public Dictionary<int, ControlPoint> ControlPointsDictionary = new Dictionary<int, ControlPoint>();

        public ModelCreator(Model model)
        {
            this.Model = model;
        }

        private bool splitToContinuousPatches =false;

        public bool SplitToContinuousPatches
        {
            get { return splitToContinuousPatches; }
            set { splitToContinuousPatches = value; }
        }

        public void CreateModelData()
        {
            if (splitToContinuousPatches)
            {
                AssignContinuousPatchesToModel();
            }else
            {
                AssignPatchesToModel();
            }            
        }

        private void AssignPatchesToModel()
        {
            int counterElementID = 0;
            int counterCPID = 0;
            for (int patchID = 0; patchID < NumberOfPatches; patchID++)
            {
                Patch patch = new Patch()
                {
                    NumberOfDimensions = this.NumberOfDimensions,
                    ID = patchID,
                    DegreeKsi = DegreeKsiDictionary[patchID],
                    DegreeHeta = DegreeHetaDictionary[patchID],
                    DegreeZeta = (NumberOfDimensions == 2) ? 0 : DegreeZetaDictionary[patchID],
                    NumberOfControlPointsKsi = NumberOfControlPointsKsiDictionary[patchID],
                    NumberOfControlPointsHeta = NumberOfControlPointsHetaDictionary[patchID],
                    NumberOfControlPointsZeta = (NumberOfDimensions == 2) ? 0 : NumberOfControlPointsZetaDictionary[patchID],
                    Material = this.Material,
                    Thickness = (NumberOfDimensions == 2) ? this.Thickness : 0,
                    KnotValueVectorKsi = new Vector(KnotValueVectorsKsiDictionary[patchID]),
                    KnotValueVectorHeta = new Vector(KnotValueVectorsHetaDictionary[patchID]),
                    KnotValueVectorZeta = (NumberOfDimensions == 2) ? null : new Vector(KnotValueVectorsZetaDictionary[patchID]),
                };
                
                for (int j = 0; j < ControlPointIDsDictionary[patchID].Length; j++)
                    patch.ControlPointsDictionary.Add(j,ControlPointsDictionary[ControlPointIDsDictionary[patchID][j]]);
                patch.CreatePatchData();
                foreach (var element in patch.ElementsDictionary.Values)
                    Model.ElementsDictionary.Add(counterElementID++, element);
                foreach (var controlPoint in patch.ControlPointsDictionary.Values)
                    Model.ControlPointsDictionary.Add(counterCPID++, controlPoint);

                this.Model.PatchesDictionary.Add(patchID, patch);
            }
            //Model.ConnectDataStructures();
        }

        private void AssignContinuousPatchesToModel()
        {
            int counterPatch = 0;
            for (int patchID = 0; patchID < NumberOfPatches; patchID++)
            {
                #region FindSubPatches
                var tupleKsi = DetectSubPatches(new Vector(KnotValueVectorsKsiDictionary[patchID]), DegreeKsiDictionary[patchID]);
                int subpatchesKsi = tupleKsi.Item1;
                Dictionary<int, Vector> subKnotVectorsKsi = tupleKsi.Item2;

                var tupleHeta = DetectSubPatches(new Vector(KnotValueVectorsHetaDictionary[patchID]), DegreeHetaDictionary[patchID]);
                int subpatchesHeta = tupleHeta.Item1;
                Dictionary<int, Vector> subKnotVectorsHeta = tupleHeta.Item2;

                int subpatchesZeta = 1;
                Tuple<int, Dictionary<int, Vector>> tupleZeta=new Tuple<int, Dictionary<int, Vector>>(1,null);
                Dictionary<int, Vector> subKnotVectorsZeta = new Dictionary<int, Vector>();
                if (this.NumberOfDimensions==3)
                {
                    tupleZeta = DetectSubPatches(new Vector(KnotValueVectorsZetaDictionary[patchID]), DegreeZetaDictionary[patchID]);
                    subpatchesZeta = tupleZeta.Item1;
                    subKnotVectorsZeta = tupleZeta.Item2;
                }
                #endregion
                var controlPointIDs = FindControlPointsOfEachSubPatch(patchID,tupleKsi, tupleHeta, tupleZeta);
        


                for (int i = 0; i < subpatchesKsi; i++)
                {
                    for (int j = 0; j < subpatchesHeta; j++)
                    {
                        for (int k = 0; k < subpatchesZeta; k++)
                        {
                            Patch patch = new Patch()
                            {
                                NumberOfDimensions = this.NumberOfDimensions,
                                ID = counterPatch,
                                DegreeKsi = DegreeKsiDictionary[patchID],
                                DegreeHeta = DegreeHetaDictionary[patchID],
                                DegreeZeta = (NumberOfDimensions == 2) ? 0 : DegreeZetaDictionary[patchID],
                                NumberOfControlPointsKsi = subKnotVectorsKsi[i].Length- DegreeKsiDictionary[patchID]-1,
                                NumberOfControlPointsHeta = subKnotVectorsHeta[j].Length - DegreeHetaDictionary[patchID] - 1,
                                NumberOfControlPointsZeta = (NumberOfDimensions == 2) ? 0 : subKnotVectorsZeta[k].Length - DegreeZetaDictionary[patchID] - 1,
                                Material = this.Material,
                                Thickness = (NumberOfDimensions == 2) ? this.Thickness : 0,
                                KnotValueVectorKsi = subKnotVectorsKsi[i],
                                KnotValueVectorHeta = subKnotVectorsHeta[j],
                                KnotValueVectorZeta = (NumberOfDimensions == 2) ? null : subKnotVectorsZeta[k]
                            };

                            for (int m = 0; m < controlPointIDs[i,j,k].Length; m++)
                            {
                                patch.ControlPointsDictionary.Add(m, ControlPointsDictionary[ControlPointIDsDictionary[patchID][controlPointIDs[i, j, k][m]]]);
                            }
                            patch.CreatePatchData();
                            this.Model.PatchesDictionary.Add(counterPatch++, patch);
                        }
                    }
                }
            }
            //Model.ConnectDataStructures();
        }

        private Tuple<int, Dictionary<int, Vector>> DetectSubPatches(Vector knotValueVector, int degree)
        {
            Dictionary<int, Vector> SubKnotVectors = new Dictionary<int, Vector>();
            Vector[] result =knotValueVector.RemoveDuplicatesFindMultiplicity();
            Vector singleValues = result[0];
            Vector multiplicity = result[1];

            int initialKnotVectorPosition = 0;
            int endingKnotVectorPosition = 0;
            int counterPatch = 0;
            for (int i = 0; i < singleValues.Length-1; i++)
            {
                if (multiplicity[i + 1] - multiplicity[i] + 1 == degree)
                {
                    endingKnotVectorPosition = i + (int)multiplicity[i];
                    if (initialKnotVectorPosition==0)
                    {
                        int length = endingKnotVectorPosition - initialKnotVectorPosition + 1+degree;
                        Vector subKnotVector = new Vector(length);
                        for (int j = 0 ; j < endingKnotVectorPosition- initialKnotVectorPosition + 1; j++)
                            subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j];
                        for (int j = endingKnotVectorPosition - initialKnotVectorPosition + 1; j < length; j++)
                            subKnotVector[j] = knotValueVector[endingKnotVectorPosition];
                        SubKnotVectors.Add(counterPatch++, subKnotVector);
                        initialKnotVectorPosition = endingKnotVectorPosition;
                    }
                    else
                    {
                        int length = endingKnotVectorPosition - initialKnotVectorPosition + 2 + degree;
                        Vector subKnotVector = new Vector(length);
                        subKnotVector[0]= knotValueVector[initialKnotVectorPosition];
                        for (int j = 1; j < endingKnotVectorPosition - initialKnotVectorPosition + 2; j++)
                            subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j-1];
                        for (int j = endingKnotVectorPosition - initialKnotVectorPosition + 2; j < length; j++)
                            subKnotVector[j] = knotValueVector[endingKnotVectorPosition];
                        SubKnotVectors.Add(counterPatch++, subKnotVector);
                        initialKnotVectorPosition = endingKnotVectorPosition;
                    }
                }
                else if (i == singleValues.Length - 2&& multiplicity[i + 1] - multiplicity[i] + 1 != degree)
                {
                    endingKnotVectorPosition = knotValueVector.Length - 1;
                    if (initialKnotVectorPosition == 0)
                    {
                        SubKnotVectors.Add(counterPatch++, knotValueVector);
                    }else
                    {
                        int length = endingKnotVectorPosition - initialKnotVectorPosition + 2;
                        Vector subKnotVector = new Vector(length);
                        subKnotVector[0] = knotValueVector[initialKnotVectorPosition];
                        for (int j = 1; j < endingKnotVectorPosition - initialKnotVectorPosition + 2; j++)
                            subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j-1];
                        SubKnotVectors.Add(counterPatch++, subKnotVector);
                    }
                }
                else if(initialKnotVectorPosition== singleValues.Length - degree)
                {
                    endingKnotVectorPosition = knotValueVector.Length - 1;
                    int length = endingKnotVectorPosition - initialKnotVectorPosition + 2;
                    Vector subKnotVector = new Vector(length);
                    subKnotVector[0] = knotValueVector[initialKnotVectorPosition];
                    for (int j = 1; j < endingKnotVectorPosition - initialKnotVectorPosition + 2; j++)
                        subKnotVector[j] = knotValueVector[initialKnotVectorPosition + j-1];
                    SubKnotVectors.Add(counterPatch++, subKnotVector);
                    break;
                }                
            }
            return new Tuple<int, Dictionary<int, Vector>>(counterPatch, SubKnotVectors);
        }

        private int[,,][] FindControlPointsOfEachSubPatch(int patchID,Tuple<int, Dictionary<int, Vector>> tupleKsi, Tuple<int, Dictionary<int, Vector>> tupleHeta, Tuple<int, Dictionary<int, Vector>> tupleZeta)
        {
            int[,,][] controlPointIDs;
            if (NumberOfDimensions == 2)
            {
                int[][] axisKsiControlPointIDs = CalculateAxisControlPointsIDs(tupleKsi, DegreeKsiDictionary, patchID);
                int[][] axisHetaControlPointIDs = CalculateAxisControlPointsIDs(tupleHeta, DegreeHetaDictionary, patchID);

                int numberOfSubpatchesKsi = axisKsiControlPointIDs.Length;
                int numberOfSubpatchesHeta = axisHetaControlPointIDs.Length;

                controlPointIDs = new int[numberOfSubpatchesKsi, numberOfSubpatchesHeta, 1][];

                for (int i = 0; i < numberOfSubpatchesKsi; i++)
                {
                    for (int j = 0; j < numberOfSubpatchesHeta; j++)
                    {
                        int numberOfSubPatchCP = (tupleKsi.Item2[i].Length - DegreeKsiDictionary[patchID] - 1) *
                            (tupleHeta.Item2[j].Length - DegreeHetaDictionary[patchID] - 1);
                        controlPointIDs[i, j, 0] = new int[numberOfSubPatchCP];
                        for (int m = 0; m < axisKsiControlPointIDs[i].Length; m++)
                        {
                            for (int n = 0; n < axisHetaControlPointIDs[j].Length; n++)
                            {
                                int indexLocal = m * axisHetaControlPointIDs[j].Length + n ;
                                int indexGlobal = axisKsiControlPointIDs[i][m] * NumberOfControlPointsHetaDictionary[patchID]  +
                                    axisHetaControlPointIDs[j][n];
                                controlPointIDs[i, j, 0][indexLocal] = indexGlobal;
                            }
                        }
                    }
                }
            }
            else
            {
                int[][] axisKsiControlPointIDs = CalculateAxisControlPointsIDs(tupleKsi,DegreeKsiDictionary,patchID);
                int[][] axisHetaControlPointIDs = CalculateAxisControlPointsIDs(tupleHeta, DegreeHetaDictionary, patchID);
                int[][] axisZetaControlPointIDs = CalculateAxisControlPointsIDs(tupleZeta, DegreeZetaDictionary, patchID);

                //int numberOfSubpatchesKsi = tupleKsi.Item1;
                //Dictionary<int, Vector> subKnotVectorsKsi = tupleKsi.Item2;
                //int[][] axisKsiControlPointIDs = new int[numberOfSubpatchesKsi][];
                //int counterCPKsi = 0;
                //for (int i = 0; i < numberOfSubpatchesKsi; i++)
                //{
                //    int length = subKnotVectorsKsi[i].Length - DegreeKsiDictionary[patchID] - 1;
                //    axisKsiControlPointIDs[i] = new int[length];
                //    for (int j = 0; j < length; j++)
                //    {
                //        axisKsiControlPointIDs[i][j] = counterCPKsi - i;
                //        counterCPKsi++;
                //    }
                //}

                //int numberOfSubpatchesHeta = tupleHeta.Item1;
                //Dictionary<int, Vector> subKnotVectorsHeta = tupleHeta.Item2;
                //int[][] axisHetaControlPointIDs = new int[numberOfSubpatchesHeta][];
                //int counterCPHeta = 0;
                //for (int i = 0; i < numberOfSubpatchesHeta; i++)
                //{
                //    int length = subKnotVectorsHeta[i].Length - DegreeHetaDictionary[patchID] - 1;
                //    axisHetaControlPointIDs[i] = new int[length];
                //    for (int j = 0; j < length; j++)
                //    {
                //        axisHetaControlPointIDs[i][j] = counterCPHeta - i;
                //        counterCPHeta++;
                //    }
                //}
                //int numberOfSubpatchesZeta = tupleZeta.Item1;
                //Dictionary<int, Vector> subKnotVectorsZeta = tupleZeta.Item2;
                //int[][] axisZetaControlPointIDs = new int[numberOfSubpatchesZeta][];
                //int counterCPZeta = 0;
                //for (int i = 0; i < numberOfSubpatchesZeta; i++)
                //{
                //    int length = subKnotVectorsZeta[i].Length - DegreeZetaDictionary[patchID] - 1;
                //    axisZetaControlPointIDs[i] = new int[length];
                //    for (int j = 0; j < length; j++)
                //    {
                //        axisZetaControlPointIDs[i][j] = counterCPZeta - i;
                //        counterCPKsi++;
                //    }
                //}
                int numberOfSubpatchesKsi = axisKsiControlPointIDs.Length;
                int numberOfSubpatchesHeta = axisHetaControlPointIDs.Length;
                int numberOfSubpatchesZeta = (NumberOfDimensions==3)? axisZetaControlPointIDs.Length:1;

                controlPointIDs = new int[numberOfSubpatchesKsi, numberOfSubpatchesHeta, numberOfSubpatchesZeta][];

                for (int i = 0; i < numberOfSubpatchesKsi; i++)
                {
                    for (int j = 0; j < numberOfSubpatchesHeta; j++)
                    {
                        for (int k = 0; k < numberOfSubpatchesZeta; k++)
                        {
                            int numberOfSubPatchCP = (tupleKsi.Item2[i].Length - DegreeKsiDictionary[patchID] - 1) *
                                (tupleHeta.Item2[j].Length - DegreeHetaDictionary[patchID] - 1) *
                                (tupleZeta.Item2[k].Length - DegreeZetaDictionary[patchID] - 1);
                            controlPointIDs[i, j, k] = new int[numberOfSubPatchCP];
                            for (int m = 0; m < axisKsiControlPointIDs[i].Length; m++)
                            {
                                for (int n = 0; n < axisHetaControlPointIDs[j].Length; n++)
                                {
                                    for (int p = 0; p < axisZetaControlPointIDs[k].Length; p++)
                                    {
                                        int indexLocal = m * axisHetaControlPointIDs[j].Length * axisZetaControlPointIDs[k].Length +
                                            n * axisZetaControlPointIDs[k].Length + p;
                                        int indexGlobal = axisKsiControlPointIDs[i][m] * NumberOfControlPointsHetaDictionary[patchID] * NumberOfControlPointsZetaDictionary[patchID] +
                                            axisHetaControlPointIDs[j][n] * NumberOfControlPointsZetaDictionary[patchID] + axisZetaControlPointIDs[k][p];
                                        controlPointIDs[i, j, k][indexLocal] = indexGlobal;
                                    }
                                }
                            }

                        }
                    }
                }

            }
            

            return controlPointIDs;

        }

        private int[][] CalculateAxisControlPointsIDs(Tuple<int, Dictionary<int, Vector>> tupleAxis,Dictionary<int,int> DegreeDictionary, int patchID)
        {
            int numberOfSubpatches = tupleAxis.Item1;
            Dictionary<int, Vector> subKnotVectorsKsi = tupleAxis.Item2;
            int[][] axisControlPointIDs = new int[numberOfSubpatches][];
            int counterCPKsi = 0;
            for (int i = 0; i < numberOfSubpatches; i++)
            {
                int length = subKnotVectorsKsi[i].Length - DegreeDictionary[patchID] - 1;
                axisControlPointIDs[i] = new int[length];
                for (int j = 0; j < length; j++)
                {
                    axisControlPointIDs[i][j] = counterCPKsi - i;
                    counterCPKsi++;
                }
            }
            return axisControlPointIDs;
        }

        //private int[,,][] CalculateSubPatchControlPointIDs(int[][] axisKsiControlPointIDs, int[][] axisHetaControlPointIDs, int[][] axisZetaControlPointIDs, int patchID)
        //{
        //    int numberOfSubpatchesKsi = axisKsiControlPointIDs.Length;
        //    int numberOfSubpatchesHeta = axisHetaControlPointIDs.Length;
        //    int numberOfSubpatchesZeta = (NumberOfDimensions==3)? axisZetaControlPointIDs.Length:1;
        //    int[,,][] controlPointIDs = new int[axisKsiControlPointIDs[0].Length, numberOfSubpatchesHeta, numberOfSubpatchesZeta][];

        //    for (int i = 0; i < numberOfSubpatchesKsi; i++)
        //    {
        //        for (int j = 0; j < numberOfSubpatchesHeta; j++)
        //        {
        //            for (int k = 0; k < numberOfSubpatchesZeta; k++)
        //            {
        //                int numberOfSubPatchCP = (subKnotVectorsKsi[i].Length - DegreeKsiDictionary[patchID] - 1) *
        //                    (subKnotVectorsHeta[j].Length - DegreeHetaDictionary[patchID] - 1) *
        //                    (subKnotVectorsZeta[k].Length - DegreeZetaDictionary[patchID] - 1);
        //                controlPointIDs[i, j, k] = new int[numberOfSubPatchCP];
        //                for (int m = 0; m < axisKsiControlPointIDs[i].Length; m++)
        //                {
        //                    for (int n = 0; n < axisHetaControlPointIDs[j].Length; n++)
        //                    {
        //                        for (int p = 0; p < axisZetaControlPointIDs[k].Length; p++)
        //                        {
        //                            int indexLocal = m * axisHetaControlPointIDs[j].Length * axisZetaControlPointIDs[k].Length +
        //                                n * axisZetaControlPointIDs[k].Length + p;
        //                            int indexGlobal = axisKsiControlPointIDs[i][m] * NumberOfControlPointsHetaDictionary[patchID] * NumberOfControlPointsZetaDictionary[patchID] +
        //                                axisHetaControlPointIDs[j][n] * NumberOfControlPointsZetaDictionary[patchID] + axisZetaControlPointIDs[k][p];
        //                            controlPointIDs[i, j, k][indexLocal] = indexGlobal;
        //                        }
        //                    }
        //                }

        //            }
        //        }
        //    }
        //}



    }
}
