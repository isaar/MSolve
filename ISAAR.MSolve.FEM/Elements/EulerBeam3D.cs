using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.FEM.Embedding;
using System.IO;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using IEmbeddedElement = ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
using ISAAR.MSolve.Discretization;

namespace ISAAR.MSolve.FEM.Elements
{
    public class EulerBeam3D : IStructuralFiniteElement, IEmbeddedElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[6] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly double youngModulus;
        private readonly double poissonRatio;
        private readonly List<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>();
        private const int hostDofsPerNode = 3;
        private const int embeddedDofsPerNode = 6;
        private const int commonDofsPerNode = 3;
        //private Matrix2D<double> transformation;
        private int noOfDOFs = 12;
        private DOFType[][] dofsWhenNoRotations = null;
        private List<Element> hostElementList;
        private bool[] isNodeEmbedded;
        private readonly Node[][] rotNodes = new Node[2][];
        private Matrix2D rotTransformation;
        private IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }
        public double SectionArea { get; set; }
        public double RayleighAlpha { get; set; }
        public double RayleighBeta { get; set; }
        //public double MomentOfInertiaX { get; set; }
        public double MomentOfInertiaY { get; set; }
        public double MomentOfInertiaZ { get; set; }
        public double MomentOfInertiaPolar { get; set; }
        public IList<EmbeddedNode> EmbeddedNodes { get { return embeddedNodes; } }

        public EulerBeam3D(double youngModulus, double poissonRatio)
        {
            this.youngModulus = youngModulus;
            this.poissonRatio = poissonRatio;
        }

        public EulerBeam3D(double youngModulus, double poissonRatio, Node[] rot1Nodes, Node[] rot2Nodes)
            : this(youngModulus, poissonRatio)
        {
            if (rot1Nodes != null && rot1Nodes.Length != 4)
                throw new ArgumentException("Dependent nodes quantity for rotation1 has to be four.");
            if (rot2Nodes != null && rot2Nodes.Length != 4)
                throw new ArgumentException("Dependent nodes quantity for rotation2 has to be four.");
            rotNodes[0] = rot1Nodes;
            rotNodes[1] = rot2Nodes;

            InitializeDOFsWhenNoRotations();
        }

        public EulerBeam3D(double youngModulus, double poissonRatio, IElementDOFEnumerator dofEnumerator) : this(youngModulus, poissonRatio)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public EulerBeam3D(double youngModulus, double poissonRatio, Node[] rot1Nodes, Node[] rot2Nodes, IElementDOFEnumerator dofEnumerator)
            : this(youngModulus, poissonRatio, rot1Nodes, rot2Nodes)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        private void InitializeDOFsWhenNoRotations()
        {
            if (rotNodes[0] == null && rotNodes[1] == null) return;

            DOFType[] translationalDOFTypes = new DOFType[3] { DOFType.X, DOFType.Y, DOFType.Z };
            dofsWhenNoRotations = new DOFType[][] { translationalDOFTypes, translationalDOFTypes,
                translationalDOFTypes, translationalDOFTypes, translationalDOFTypes, translationalDOFTypes,
                translationalDOFTypes, translationalDOFTypes, translationalDOFTypes, translationalDOFTypes };
            noOfDOFs = 30;

            if (rotNodes[0] == null)
            {
                dofsWhenNoRotations = new DOFType[][] { nodalDOFTypes, translationalDOFTypes, translationalDOFTypes,
                translationalDOFTypes, translationalDOFTypes, translationalDOFTypes };
                noOfDOFs = 21;
            }

            if (rotNodes[1] == null)
            {
                dofsWhenNoRotations = new DOFType[][] { translationalDOFTypes, nodalDOFTypes, translationalDOFTypes,
                translationalDOFTypes, translationalDOFTypes, translationalDOFTypes };
                noOfDOFs = 21;
            }
        }

        private void CalculateRotTranformation(IElement element)
        {
            if (rotNodes[0] == null && rotNodes[1] == null)
            {
                rotTransformation = new Matrix2D(12, 12);
                for (int i = 0; i < 12; i++) rotTransformation[i, i] = 1;
                return;
            }

            int[] transMatrixRows = new int[] { 3, 9 };
            int[] transMatrixCols = new int[] { 6, 18 };
            int[] transMatrixColRows = new int[] { 0, 3 };
            int nonRotationalDOFs = 24;
            //            int nonRotationalDOFs = 0;
            if (rotNodes[0] == null)
            {
                nonRotationalDOFs = 13;
                transMatrixRows = new int[] { -1, 9 };
                transMatrixCols = new int[] { -1, 9 };
                transMatrixColRows = new int[] { -1, 6 };
            }
            if (rotNodes[1] == null)
            {
                nonRotationalDOFs = 13;
                transMatrixRows = new int[] { 3, -1 };
                transMatrixCols = new int[] { 9, -1 };
                transMatrixColRows = new int[] { 0, -1 };
            }

            for (int i = 0; i < 2; i++)
                if (rotNodes[i] == null)
                    nonRotationalDOFs += 2;
            //else
            //    nonRotationalDOFs += 12;

            rotTransformation = new Matrix2D(12, nonRotationalDOFs + 6);
            rotTransformation[0, 0] = 1;
            rotTransformation[1, 1] = 1;
            rotTransformation[2, 2] = 1;
            if (rotNodes[0] == null)
            {
                rotTransformation[3, 3] = 1;
                rotTransformation[4, 4] = 1;
                rotTransformation[5, 5] = 1;
                rotTransformation[6, 6] = 1;
                rotTransformation[7, 7] = 1;
                rotTransformation[8, 8] = 1;
            }
            else if (rotNodes[1] == null)
            {
                rotTransformation[6, 3] = 1;
                rotTransformation[7, 4] = 1;
                rotTransformation[8, 5] = 1;
                rotTransformation[9, 6] = 1;
                rotTransformation[10, 7] = 1;
                rotTransformation[11, 8] = 1;
            }
            else
            {
                rotTransformation[6, 3] = 1;
                rotTransformation[7, 4] = 1;
                rotTransformation[8, 5] = 1;
            }

            double[][] rotDifsX = new double[2][];
            double[][] rotDifsY = new double[2][];
            double[][] rotDifsZ = new double[2][];
            double[][] lengthsSquared = new double[2][];
            for (int i = 0; i < 2; i++)
            {
                if (rotNodes[i] == null) continue;

                rotDifsX[i] = new double[]
                {
                    rotNodes[i][0].X - element.INodes[i].X,
                    rotNodes[i][1].X - element.INodes[i].X,
                    rotNodes[i][2].X - element.INodes[i].X,
                    rotNodes[i][3].X - element.INodes[i].X
                };
                rotDifsY[i] = new double[]
                {
                    rotNodes[i][0].Y - element.INodes[i].Y,
                    rotNodes[i][1].Y - element.INodes[i].Y,
                    rotNodes[i][2].Y - element.INodes[i].Y,
                    rotNodes[i][3].Y - element.INodes[i].Y
                };
                rotDifsZ[i] = new double[]
                {
                    rotNodes[i][0].Z - element.INodes[i].Z,
                    rotNodes[i][1].Z - element.INodes[i].Z,
                    rotNodes[i][2].Z - element.INodes[i].Z,
                    rotNodes[i][3].Z - element.INodes[i].Z
                };

                lengthsSquared[i] = new double[]
                {
                    Math.Pow(rotDifsX[i][0], 2) + Math.Pow(rotDifsY[i][0], 2) + Math.Pow(rotDifsZ[i][0], 2),
                    Math.Pow(rotDifsX[i][1], 2) + Math.Pow(rotDifsY[i][1], 2) + Math.Pow(rotDifsZ[i][1], 2),
                    Math.Pow(rotDifsX[i][2], 2) + Math.Pow(rotDifsY[i][2], 2) + Math.Pow(rotDifsZ[i][2], 2),
                    Math.Pow(rotDifsX[i][3], 2) + Math.Pow(rotDifsY[i][3], 2) + Math.Pow(rotDifsZ[i][3], 2)
                };
            }

            for (int i = 0; i < 2; i++)
            {
                if (rotNodes[i] == null) continue;

                int r = transMatrixRows[i];
                int c = transMatrixCols[i];
                int cr = transMatrixColRows[i];

                rotTransformation[r + 0, cr + 1] = rotDifsY[i][3] / lengthsSquared[i][3] -
                    rotDifsY[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 0, cr + 2] = rotDifsZ[i][3] / lengthsSquared[i][3] -
                    rotDifsZ[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 0, c + 4] = rotDifsY[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 0, c + 5] = rotDifsZ[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 0, c + 10] = -rotDifsY[i][3] / lengthsSquared[i][3];
                rotTransformation[r + 0, c + 11] = -rotDifsZ[i][3] / lengthsSquared[i][3];

                rotTransformation[r + 1, cr + 0] = rotDifsX[i][2] / lengthsSquared[i][2] -
                    rotDifsX[i][0] / lengthsSquared[i][0] + rotDifsZ[i][3] / lengthsSquared[i][3] -
                    rotDifsZ[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 1, cr + 2] = rotDifsZ[i][2] / lengthsSquared[i][2] -
                    rotDifsZ[i][0] / lengthsSquared[i][0] + rotDifsX[i][3] / lengthsSquared[i][3] -
                    rotDifsX[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 1, c + 0] = rotDifsX[i][0] / lengthsSquared[i][0];
                rotTransformation[r + 1, c + 2] = rotDifsZ[i][0] / lengthsSquared[i][0];
                rotTransformation[r + 1, c + 3] = rotDifsX[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 1, c + 5] = rotDifsZ[i][1] / lengthsSquared[i][1];
                rotTransformation[r + 1, c + 6] = -rotDifsX[i][2] / lengthsSquared[i][2];
                rotTransformation[r + 1, c + 8] = -rotDifsZ[i][2] / lengthsSquared[i][2];
                rotTransformation[r + 1, c + 9] = -rotDifsX[i][3] / lengthsSquared[i][3];
                rotTransformation[r + 1, c + 11] = -rotDifsZ[i][3] / lengthsSquared[i][3];

                rotTransformation[r + 2, cr + 0] = rotDifsX[i][2] / lengthsSquared[i][2] -
                    rotDifsX[i][0] / lengthsSquared[i][0];
                rotTransformation[r + 2, cr + 1] = rotDifsY[i][2] / lengthsSquared[i][2] -
                    rotDifsY[i][0] / lengthsSquared[i][0];
                rotTransformation[r + 2, c + 0] = rotDifsX[i][0] / lengthsSquared[i][0];
                rotTransformation[r + 2, c + 1] = rotDifsY[i][0] / lengthsSquared[i][0];
                rotTransformation[r + 2, c + 6] = -rotDifsX[i][2] / lengthsSquared[i][2];
                rotTransformation[r + 2, c + 7] = -rotDifsY[i][2] / lengthsSquared[i][2];

                //rotTransformation[r + 0, cr + 1] = rotDifsY[i][3] / lengthsSquared[i][3] +
                //    rotDifsY[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 0, cr + 2] = rotDifsZ[i][3] / lengthsSquared[i][3] +
                //    rotDifsZ[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 0, c + 4] = rotDifsY[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 0, c + 5] = rotDifsZ[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 0, c + 10] = rotDifsY[i][3] / lengthsSquared[i][3];
                //rotTransformation[r + 0, c + 11] = rotDifsZ[i][3] / lengthsSquared[i][3];

                //rotTransformation[r + 1, cr + 0] = rotDifsX[i][2] / lengthsSquared[i][2] +
                //    rotDifsX[i][0] / lengthsSquared[i][0] + rotDifsZ[i][3] / lengthsSquared[i][3] +
                //    rotDifsZ[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 1, cr + 2] = rotDifsZ[i][2] / lengthsSquared[i][2] +
                //    rotDifsZ[i][0] / lengthsSquared[i][0] + rotDifsX[i][3] / lengthsSquared[i][3] +
                //    rotDifsX[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 1, c + 0] = rotDifsX[i][0] / lengthsSquared[i][0];
                //rotTransformation[r + 1, c + 2] = rotDifsZ[i][0] / lengthsSquared[i][0];
                //rotTransformation[r + 1, c + 3] = rotDifsX[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 1, c + 5] = rotDifsZ[i][1] / lengthsSquared[i][1];
                //rotTransformation[r + 1, c + 6] = rotDifsX[i][2] / lengthsSquared[i][2];
                //rotTransformation[r + 1, c + 8] = rotDifsZ[i][2] / lengthsSquared[i][2];
                //rotTransformation[r + 1, c + 9] = rotDifsX[i][3] / lengthsSquared[i][3];
                //rotTransformation[r + 1, c + 11] = rotDifsZ[i][3] / lengthsSquared[i][3];

                //rotTransformation[r + 2, cr + 0] = rotDifsX[i][2] / lengthsSquared[i][2] +
                //    rotDifsX[i][0] / lengthsSquared[i][0];
                //rotTransformation[r + 2, cr + 1] = rotDifsY[i][2] / lengthsSquared[i][2] +
                //    rotDifsY[i][0] / lengthsSquared[i][0];
                //rotTransformation[r + 2, c + 0] = rotDifsX[i][0] / lengthsSquared[i][0];
                //rotTransformation[r + 2, c + 1] = rotDifsY[i][0] / lengthsSquared[i][0];
                //rotTransformation[r + 2, c + 6] = rotDifsX[i][2] / lengthsSquared[i][2];
                //rotTransformation[r + 2, c + 7] = rotDifsY[i][2] / lengthsSquared[i][2];
            }
        }

        #region IElementType Members

        public int ID
        {
            get { return 2; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        private IList<Tuple<Node, IList<DOFType>>> GetDOFTypesInternal(Element element)
        {
            if (element == null) throw new ArgumentException();

            var hostDOFTypes = new List<DOFType>();
            foreach (var node in element.Nodes)
            {
                var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
                if (embeddedNode != null)
                    hostDOFTypes.AddRange(embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(null).SelectMany(x => x));
            }
            hostDOFTypes = hostDOFTypes.Distinct().ToList();

            var d = new Dictionary<Node, IList<DOFType>>();
            var l = new List<Tuple<Node, IList<DOFType>>>();
            foreach (var node in element.Nodes)
            {
                var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
                //if (node.EmbeddedInElement == null)
                if (embeddedNode == null)
                {
                    var nodeDofs = new List<DOFType>();
                    nodeDofs.AddRange(nodalDOFTypes.Except(hostDOFTypes));
                    d.Add(node, nodeDofs);
                    l.Add(new Tuple<Node, IList<DOFType>>(node, nodeDofs));
                }
                else
                {
                    //d.AddRange(node.EmbeddedInElement.ElementType.GetDOFTypes(null));
                    var hostDOFsPerNode = embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(null);
                    for (int i = 0; i < hostDOFsPerNode.Count; i++)
                    {
                        if (!d.ContainsKey(embeddedNode.EmbeddedInElement.Nodes[i]))
                            d.Add(embeddedNode.EmbeddedInElement.Nodes[i], hostDOFsPerNode[i]);
                        l.Add(new Tuple<Node, IList<DOFType>>(embeddedNode.EmbeddedInElement.Nodes[i], hostDOFsPerNode[i]));
                    }
                }
            }

            //var hostDOFTypes = d.SelectMany(x => x.Value).Distinct();
            var uniqueDOFTypes = nodalDOFTypes.Except(hostDOFTypes).Union(hostDOFTypes.Except(nodalDOFTypes)).ToArray();
            if (uniqueDOFTypes.Length > 0)
                foreach (var node in element.Nodes)
                {
                    //if (embeddedNodes.Where(x => x.Node == node).FirstOrDefault() != null)
                    //{
                    //d.Add(node, uniqueDOFTypes);
                    //l.Add(new Tuple<Node, IList<DOFType>>(node, uniqueDOFTypes));
                    //}
                    if (!d.ContainsKey(node))
                        d.Add(node, uniqueDOFTypes);
                    else
                        d[node] = d[node].Concat(uniqueDOFTypes).ToArray();
                    l.Add(new Tuple<Node, IList<DOFType>>(node, uniqueDOFTypes));
                }

            return l;
        }

        //public IList<IList<DOFType>> GetDOFTypes(Element element)
        //{
        //    if (element == null) return dofs;

        //    var d = new List<IList<DOFType>>();
        //    var dofTypeDictionary = GetDOFTypesInternal(element);
        //    foreach (var dofTypes in dofTypeDictionary)
        //        d.Add(dofTypes.Item2);

        //    return d;
        //}

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            if (dofsWhenNoRotations == null) return dofs;
            return dofsWhenNoRotations;

            //if (element == null) return dofs;

            //var d = new List<IList<DOFType>>();
            //foreach (var node in element.Nodes)
            //{
            //    var nodeDofs = new List<DOFType>();
            //    nodeDofs.AddRange(nodalDOFTypes);
            //    d.Add(nodeDofs);
            //}
            //return d;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            var nodes = new List<Node>();
            var dofTypeDictionary = GetDOFTypesInternal(element);
            foreach (var dofType in dofTypeDictionary)
                nodes.Add(dofType.Item1);
            //foreach (var dofType in dofTypeDictionary)
            //    nodes.Add(dofType.Key);

            //foreach (var node in element.Nodes)
            //{
            //    var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
            //    //if (node.EmbeddedInElement == null)
            //    if (embeddedNode == null)
            //        nodes.Add(node);
            //    else
            //        //nodes.AddRange(node.EmbeddedInElement.Nodes);
            //        nodes.AddRange(embeddedNode.EmbeddedInElement.Nodes);
            //}

            return nodes;
        }

        private IMatrix2D StiffnessMatrixPure(IElement element)
        {
            double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            double z2 = Math.Pow(element.INodes[1].Z - element.INodes[0].Z, 2);
            double L = 1 / Math.Sqrt(x2 + y2 + z2);
            double L2 = L * L;
            double L3 = L2 * L;
            //double EIx = m.YoungModulus * MomentOfInertiaX;
            double EIy = this.youngModulus * MomentOfInertiaY;
            double EIz = this.youngModulus * MomentOfInertiaZ;
            double GJL = this.youngModulus * L * MomentOfInertiaPolar / (2 * (1 + this.poissonRatio));
            double EAL = this.youngModulus * SectionArea * L;
            var stiffnessMatrix = new SymmetricMatrix2D(new double[] { EAL, 0, 0, 0, 0, 0, -EAL, 0, 0, 0, 0, 0,
                12*EIz*L3, 0, 0, 0, 6*EIz*L2, 0, -12*EIz*L3, 0, 0, 0, 6*EIz*L2,
                12*EIy*L3, 0, -6*EIy*L2, 0, 0, 0, -12*EIy*L3, 0, -6*EIy*L2, 0,
                GJL, 0, 0, 0, 0, 0, -GJL, 0, 0,
                4*EIy*L, 0, 0, 0, 6*EIy*L2, 0, 2*EIy*L, 0,
                4*EIz*L, 0, -6*EIz*L2, 0, 0, 0, 2*EIz*L,
                EAL, 0, 0, 0, 0, 0,
                12*EIz*L3, 0, 0, 0, -6*EIz*L2,
                12*EIy*L3, 0, 6*EIy*L2, 0,
                GJL, 0, 0,
                4*EIy*L, 0,
                4*EIz*L
            });

            var refx = new double[] { 1, 1, 1 };
            var beamTransformation = new Matrix2D(12, 12);
            beamTransformation[0, 0] = (element.INodes[1].X - element.INodes[0].X) * L;
            beamTransformation[0, 1] = (element.INodes[1].Y - element.INodes[0].Y) * L;
            beamTransformation[0, 2] = (element.INodes[1].Z - element.INodes[0].Z) * L;

            //beamTransformation[2, 0] = refx[0];
            //beamTransformation[2, 1] = refx[1];
            //beamTransformation[2, 2] = refx[2];

            //beamTransformation[1, 0] = beamTransformation[2, 1] * beamTransformation[0, 2] - beamTransformation[2, 2] * beamTransformation[0, 1];
            //beamTransformation[1, 1] = beamTransformation[2, 2] * beamTransformation[0, 0] - beamTransformation[2, 0] * beamTransformation[0, 2];
            //beamTransformation[1, 2] = beamTransformation[2, 0] * beamTransformation[0, 1] - beamTransformation[2, 1] * beamTransformation[0, 0];
            beamTransformation[1, 0] = refx[1] * beamTransformation[0, 2] - refx[2] * beamTransformation[0, 1];
            beamTransformation[1, 1] = refx[2] * beamTransformation[0, 0] - refx[0] * beamTransformation[0, 2];
            beamTransformation[1, 2] = refx[0] * beamTransformation[0, 1] - refx[1] * beamTransformation[0, 0];
            double dn = 1.0 / Math.Sqrt(beamTransformation[1, 0] * beamTransformation[1, 0] + beamTransformation[1, 1] * beamTransformation[1, 1] + beamTransformation[1, 2] * beamTransformation[1, 2]);
            beamTransformation[1, 0] = beamTransformation[1, 0] * dn;
            beamTransformation[1, 1] = beamTransformation[1, 1] * dn;
            beamTransformation[1, 2] = beamTransformation[1, 2] * dn;
            beamTransformation[2, 0] = beamTransformation[0, 1] * beamTransformation[1, 2] - beamTransformation[0, 2] * beamTransformation[1, 1];
            beamTransformation[2, 1] = beamTransformation[0, 2] * beamTransformation[1, 0] - beamTransformation[0, 0] * beamTransformation[1, 2];
            beamTransformation[2, 2] = beamTransformation[0, 0] * beamTransformation[1, 1] - beamTransformation[0, 1] * beamTransformation[1, 0];

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    beamTransformation[i + 3, j + 3] = beamTransformation[i, j];
                    beamTransformation[i + 6, j + 6] = beamTransformation[i, j];
                    beamTransformation[i + 9, j + 9] = beamTransformation[i, j];
                }

            return new SymmetricMatrix2D(beamTransformation.Transpose() * stiffnessMatrix.ToMatrix2D() * beamTransformation);

            ////if (element.Nodes.Count(n => n.EmbeddedInElement != null) == 0) return stiffnessMatrix;
            //stiffnessMatrix = new SymmetricMatrix2D<double>(beamTransformation.Transpose() * stiffnessMatrix.ToMatrix2D() * beamTransformation);
            //if (embeddedNodes.Count == 0) return stiffnessMatrix;

            ////var hostElements = element.Nodes.Select(x => x.EmbeddedInElement).Distinct();
            //var size = GetElementDOFTypes(element).SelectMany(x => x).Count();
            //transformation = new Matrix2D<double>(dofs.SelectMany(d => d).Count(), size);
            //isNodeEmbedded = new bool[element.Nodes.Count];

            ////TODO : SEPARATE FROM ELEMENT!!
            ////TODO: Must match DOFs of host with embedded element
            //int row = 0;
            //int col = 0;
            //hostElementList = new List<Element>();
            //for (int i = 0; i < element.Nodes.Count; i++)
            //{
            //    var node = element.Nodes[i];
            //    var embeddedNode = embeddedNodes.Where(x => x.Node == node).FirstOrDefault();
            //    //var hostElement = node.EmbeddedInElement;
            //    Element hostElement = embeddedNode == null ? null : embeddedNode.EmbeddedInElement;
            //    if (hostElement == null)
            //    {
            //        isNodeEmbedded[i] = false;
            //        for (int j = 0; j < dofs[i].Length; j++)
            //        {
            //            transformation[row, col] = 1;
            //            row++;
            //            col++;
            //        }
            //    }
            //    else
            //    {
            //        isNodeEmbedded[i] = true;
            //        //double[] hostShapeFunctions = ((IEmbeddedHostElement)hostElement.ElementType).GetShapeFunctionsForNode(hostElement, node);
            //        double[] hostShapeFunctions = ((IEmbeddedHostElement)hostElement.ElementType).GetShapeFunctionsForNode(hostElement, embeddedNode);

            //        if (hostElementList.IndexOf(hostElement) < 0)
            //            hostElementList.Add(hostElement);
            //        else
            //            col -= hostShapeFunctions.Length * hostDofsPerNode;

            //        for (int j = 0; j < commonDofsPerNode; j++)
            //        {
            //            for (int k = 0; k < hostShapeFunctions.Length; k++)
            //                transformation[row, hostDofsPerNode * k + col + j] = hostShapeFunctions[k];
            //            row++;
            //        }
            //        row += embeddedDofsPerNode - commonDofsPerNode;
            //        col += hostShapeFunctions.Length * hostDofsPerNode;
            //    }
            //}

            //// Add identity matrix
            //int index = 0;
            //if (isNodeEmbedded[0])
            //{
            //    transformation[3, col] = 1;
            //    transformation[4, col + 1] = 1;
            //    transformation[5, col + 2] = 1;
            //    index += 3;
            //}
            //if (isNodeEmbedded[1])
            //{
            //    transformation[9, col + index] = 1;
            //    transformation[10, col + index + 1] = 1;
            //    transformation[11, col + index + 2] = 1;
            //}

            //var transformedMatrix = new SymmetricMatrix2D<double>(transformation.Transpose() * stiffnessMatrix.ToMatrix2D() * transformation);
            ////var sw = File.CreateText(@"d:\BeamTransformed.txt");
            ////for (int i = 0; i < 54; i++)
            ////{
            ////    var s = string.Empty;
            ////    for (int j = 0; j < 54; j++)
            ////        s += transformedMatrix[i, j].ToString() + ";";
            ////    sw.WriteLine(s);
            ////}
            ////sw.Close();
            //return transformedMatrix;
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            CalculateRotTranformation(element);
            return dofEnumerator.GetTransformedMatrix(new SymmetricMatrix2D(rotTransformation.Transpose() * ((SymmetricMatrix2D)StiffnessMatrixPure(element)).ToMatrix2D() * rotTransformation));
        }

        public IMatrix2D MassMatrix(IElement element)
        {
            double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            double z2 = Math.Pow(element.INodes[1].Z - element.INodes[0].Z, 2);
            double L = 1d / Math.Sqrt(x2 + y2 + z2);
            //double halfMass = 0.5 * Density * SectionArea * L;

            //var massMatrix = new SymmetricMatrix2D<double>(new double[] { halfMass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //    halfMass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //    halfMass, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //    halfMass, 0, 0, 0, 0, 0, 0, 0, 0,
            //    halfMass, 0, 0, 0, 0, 0, 0, 0,
            //    halfMass, 0, 0, 0, 0, 0, 0,
            //    halfMass, 0, 0, 0, 0, 0,
            //    halfMass, 0, 0, 0, 0, 
            //    halfMass, 0, 0, 0,
            //    halfMass, 0, 0,
            //    halfMass, 0,
            //    halfMass
            //});
            double halfMass = Density * SectionArea / L / 6d;
            var massMatrix = new SymmetricMatrix2D(new double[] { halfMass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                halfMass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                halfMass, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0,
                halfMass, 0, 0, 0, 0, 0,
                halfMass, 0, 0, 0, 0,
                halfMass, 0, 0, 0,
                0, 0, 0,
                0, 0,
                0
            });

            var refx = new double[] { 1, 1, 1 };
            var beamTransformation = new Matrix2D(12, 12);
            beamTransformation[0, 0] = (element.INodes[1].X - element.INodes[0].X) * L;
            beamTransformation[0, 1] = (element.INodes[1].Y - element.INodes[0].Y) * L;
            beamTransformation[0, 2] = (element.INodes[1].Z - element.INodes[0].Z) * L;

            beamTransformation[1, 0] = refx[1] * beamTransformation[0, 2] - refx[2] * beamTransformation[0, 1];
            beamTransformation[1, 1] = refx[2] * beamTransformation[0, 0] - refx[0] * beamTransformation[0, 2];
            beamTransformation[1, 2] = refx[0] * beamTransformation[0, 1] - refx[1] * beamTransformation[0, 0];
            double dn = 1.0 / Math.Sqrt(beamTransformation[1, 0] * beamTransformation[1, 0] + beamTransformation[1, 1] * beamTransformation[1, 1] + beamTransformation[1, 2] * beamTransformation[1, 2]);
            beamTransformation[1, 0] = beamTransformation[1, 0] * dn;
            beamTransformation[1, 1] = beamTransformation[1, 1] * dn;
            beamTransformation[1, 2] = beamTransformation[1, 2] * dn;
            beamTransformation[2, 0] = beamTransformation[0, 1] * beamTransformation[1, 2] - beamTransformation[0, 2] * beamTransformation[1, 1];
            beamTransformation[2, 1] = beamTransformation[0, 2] * beamTransformation[1, 0] - beamTransformation[0, 0] * beamTransformation[1, 2];
            beamTransformation[2, 2] = beamTransformation[0, 0] * beamTransformation[1, 1] - beamTransformation[0, 1] * beamTransformation[1, 0];

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                    beamTransformation[i + 3, j + 3] = beamTransformation[i, j];
                    beamTransformation[i + 6, j + 6] = beamTransformation[i, j];
                    beamTransformation[i + 9, j + 9] = beamTransformation[i, j];
                }
            CalculateRotTranformation(element);

            return dofEnumerator.GetTransformedMatrix(new SymmetricMatrix2D(rotTransformation.Transpose() * beamTransformation.Transpose() * massMatrix.ToMatrix2D() * beamTransformation * rotTransformation));
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            var m = MassMatrix(element);
            var lc = m as ILinearlyCombinable;
            lc.LinearCombination(new double[] { RayleighAlpha, RayleighBeta }, new IMatrix2D[] { MassMatrix(element), StiffnessMatrix(element) });
            return m;
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            return new Tuple<double[], double[]>(new double[6], new double[6]);
            //throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            CalculateRotTranformation(element);
            IMatrix2D stiffnessMatrix = StiffnessMatrixPure(element);
            var disps = rotTransformation * new Vector(localDisplacements);
            double[] forces = new double[disps.Length];
            stiffnessMatrix.Multiply(disps, forces);
            return forces;
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            IMatrix2D stiffnessMatrix = StiffnessMatrix(element);
            Vector disps = new Vector(localDisplacements.Length);
            double[] forces = new double[localDisplacements.Length];
            for (int i = 0; i < localDisplacements.Length; i++)
                //disps[i] = localDisplacements[i] + localdDisplacements[i];
                disps[i] = localDisplacements[i];
            stiffnessMatrix.Multiply(disps, forces);
            return forces;
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector accelerations = new Vector(noOfDOFs);
            IMatrix2D massMatrix = MassMatrix(element);

            foreach (MassAccelerationLoad load in loads)
            {
                int index = 0;
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }
            }
            double[] forces = new double[accelerations.Length];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void ClearMaterialState()
        {
        }

        public void SaveMaterialState()
        {
            //throw new NotImplementedException();
        }

        public bool MaterialModified
        {
            get { return false; }
        }

        public void ResetMaterialModified()
        {
        }

        public void ClearMaterialStresses()
        {
            //throw new NotImplementedException();
        }

        #endregion

        #region IStructuralFiniteElement Members

        //public IFiniteElementMaterial Material
        //{
        //    get { return material; }
        //}

        #endregion

        #region IEmbeddedElement Members

        public Dictionary<DOFType, int> GetInternalNodalDOFs(Element element, Node node)
        {
            int index = 0;
            foreach (var elementNode in element.Nodes)
            {
                if (node.ID == elementNode.ID)
                    break;
                index++;
            }
            if (index >= 2)
                throw new ArgumentException(String.Format("GetInternalNodalDOFs: Node {0} not found in element {1}.", node.ID, element.ID));

            return index == 0 ? new Dictionary<DOFType, int>() {
                { DOFType.X, 0 }, { DOFType.Y, 1 }, { DOFType.Z, 2 }, { DOFType.RotX, 3 }, { DOFType.RotY, 4 }, { DOFType.RotZ, 5 } } :
                new Dictionary<DOFType, int>() {
                { DOFType.X, 6 }, { DOFType.Y, 7 }, { DOFType.Z, 8 }, { DOFType.RotX, 9 }, { DOFType.RotY, 10 }, { DOFType.RotZ, 11 } };
        }

        public double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues)
        {
            //if (transformation == null)
            //    throw new InvalidOperationException("Requested embedded node values for element that has no embedded nodes.");
            //if (hostElementList == null)
            //    throw new InvalidOperationException("Requested host element list for element that has no embedded nodes.");
            //int index = hostElementList.IndexOf(hostElement);
            //if (index < 0)
            //    throw new ArgumentException("Requested host element is not inside host element list.");

            //double[] values = new double[transformation.Columns];
            //int multiplier = hostElement.ElementType.DOFEnumerator.GetDOFTypes(hostElement).SelectMany(d => d).Count();
            //int vectorIndex = 0;
            //for (int i = 0; i < index; i++)
            //    vectorIndex += isNodeEmbedded[i] ? 3 : multiplier;
            //Array.Copy(hostDOFValues, 0, values, vectorIndex, multiplier);

            //return (transformation * new Vector<double>(values)).Data;

            return dofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
        }

        #endregion

    }
}
