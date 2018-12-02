using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities.BoundaryConditions;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    public class Patch: ISubdomain
	{
        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
        private readonly Dictionary<int, ControlPoint> controlPointsDictionary = new Dictionary<int, ControlPoint>();
        private readonly Dictionary<int, Dictionary<DOFType, int>> controlPointDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private readonly Dictionary<int, Dictionary<DOFType, int>> globalControlPointsDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();
        private readonly Dictionary<int, Edge> edgesDictionary = new Dictionary<int, Edge>();
        private readonly Dictionary<int, Face> facesDictionary = new Dictionary<int, Face>();

        private double[] forces;

        #region PatchData
        public int NumberOfDimensions { get; set; }
        public double Thickness { get; set; }
        public IFiniteElementMaterial Material { get; set; }

        public int NumberOfControlPointsKsi { get; set; }
        public int NumberOfControlPointsHeta { get; set; }
        public int NumberOfControlPointsZeta { get; set; }

        public int DegreeKsi { get; set; }
        public int DegreeHeta { get; set; }
        public int DegreeZeta { get; set; }

        public IVector KnotValueVectorKsi { get; set; }
        public IVector KnotValueVectorHeta { get; set; }
        public IVector KnotValueVectorZeta { get; set; }
        #endregion

        #region Properties

        public Dictionary<int, Element> ElementsDictionary
        {
            get { return elementsDictionary; }
        }

        public Dictionary<int, ControlPoint> ControlPointsDictionary
        {
            get { return controlPointsDictionary; }
        }

        public IList<ControlPoint> ControlPoints
        {
            get { return controlPointsDictionary.Values.ToList<ControlPoint>(); }
        }

        public Dictionary<int, Dictionary<DOFType,int>> ControlPointDOFsDictionary
        {
            get { return controlPointDOFsDictionary; }
        }

        public Dictionary<int, Dictionary<DOFType, int>> GlobalControlPointDOFsDictionary
        {
            get { return globalControlPointsDOFsDictionary; }
        }

		public Dictionary<int, Dictionary<DOFType, int>> NodalDOFsDictionary
		{
			get { return controlPointDOFsDictionary; }
		}

		public Dictionary<int, Dictionary<DOFType, int>> GlobalNodalDOFsDictionary
		{
			get { return globalControlPointsDOFsDictionary; }
		}

		public double[] Forces
        {
            get { return forces; }
        }

        public bool MaterialsModified { get; set; }

        public int ID { get; set; }
        public int TotalDOFs { get; set; }

        public Dictionary<int, Edge> EdgesDictionary
        {
            get { return edgesDictionary; }
        }

        public Dictionary<int,Face> FacesDictionary
        {
            get { return facesDictionary; }
        }

		public Dictionary<int, IElement> ΙElementsDictionary
		{
			get
			{
				var a = new Dictionary<int, IElement>();
				foreach (var element in elementsDictionary.Values)
					a.Add(element.ID, element);
				return a;
			}
		}

		#endregion

		#region Data Interconnection routines

		public void EnumerateDOFs()
        {
            TotalDOFs = 0;
            Dictionary<int, List<DOFType>> controlPointDOFTypesDictionary = new Dictionary<int, List<DOFType>>();
            foreach (Element element in elementsDictionary.Values)
            {
                for (int i=0; i<element.ControlPoints.Count; i++)
                {
                    if (!controlPointDOFTypesDictionary.ContainsKey(element.ControlPoints[i].ID))
                        controlPointDOFTypesDictionary.Add(element.ControlPoints[i].ID, new List<DOFType>());
                    controlPointDOFTypesDictionary[element.ControlPoints[i].ID].AddRange(element.ElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element)[i]);
                }
            }

            foreach (ControlPoint controlPoint in controlPointsDictionary.Values)
            {
                Dictionary<DOFType, int> dofsDictionary = new Dictionary<DOFType, int>();
                foreach (DOFType dofType in controlPointDOFTypesDictionary[controlPoint.ID].Distinct<DOFType>())
                {
                    int dofID = 0;
                    foreach (DOFType constraint in controlPoint.Constrains)
                    {
                        if (constraint == dofType)
                        {
                            dofID = -1;
                            break;
                        }
                    }

                    if (dofID == 0)
                    {
                        dofID = TotalDOFs;
                        TotalDOFs++;
                    }
                    dofsDictionary.Add(dofType, dofID);
                }
                controlPointDOFsDictionary.Add(controlPoint.ID, dofsDictionary);                
            }
            forces = new double[TotalDOFs];
        }

        public void AssignGlobalControlPointDOFsFromModel(Dictionary<int, Dictionary<DOFType, int>> glodalDOFsDictionary)
        {
            foreach (int controlPointID in controlPointDOFsDictionary.Keys)
            {
                Dictionary<DOFType, int> dofTypes = controlPointDOFsDictionary[controlPointID];
                Dictionary<DOFType, int> globalDOFTypes = new Dictionary<DOFType, int>(dofTypes.Count);
                foreach (DOFType dofType in dofTypes.Keys)
                    globalDOFTypes.Add(dofType, glodalDOFsDictionary[controlPointID][dofType]);
                globalControlPointsDOFsDictionary.Add(controlPointID, globalDOFTypes);
            }
        }

        public void BuildFacesDictionary()
        {
            if (this.NumberOfDimensions<=2)
            {
                FacesDictionary.Clear();
            }else
            {
                #region FaceRight
                Face faceRight = new Face();
                faceRight.Degrees[0] = this.DegreeHeta;
                faceRight.Degrees[1] = this.DegreeZeta;
                faceRight.KnotValueVectors.Add(0,this.KnotValueVectorHeta);
                faceRight.KnotValueVectors.Add(1, this.KnotValueVectorZeta );
                faceRight.Patch = this;
                int counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
                    {
                        faceRight.ControlPointsDictionary.Add(counter++,
                            this.controlPointsDictionary[j + this.NumberOfControlPointsZeta * i]);
                    }
                }
                FacesDictionary.Add(0, faceRight);
                #endregion
                #region FaceLeft
                Face faceLeft = new Face();
                faceLeft.Degrees[0] = this.DegreeHeta;
                faceLeft.Degrees[1] = this.DegreeZeta;
                faceLeft.KnotValueVectors.Add(0, this.KnotValueVectorHeta);
                faceLeft.KnotValueVectors.Add(1, this.KnotValueVectorZeta);
                faceLeft.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
                    {
                        faceLeft.ControlPointsDictionary.Add(counter++, 
                            this.controlPointsDictionary[j + this.NumberOfControlPointsZeta * i+
                            this.NumberOfControlPointsHeta*this.NumberOfControlPointsZeta*(this.NumberOfControlPointsKsi-1)]);
                    }
                }
                FacesDictionary.Add(1, faceLeft);
                #endregion
                #region FaceBottom
                Face faceBottom = new Face();
                faceBottom.Degrees[0] = this.DegreeKsi;
                faceBottom.Degrees[1] = this.DegreeHeta;
                faceBottom.KnotValueVectors.Add(0, this.KnotValueVectorKsi);
                faceBottom.KnotValueVectors.Add(1, this.KnotValueVectorHeta);
                faceBottom.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    for (int j = 0; j < this.NumberOfControlPointsHeta; j++)
                    {
                        faceBottom.ControlPointsDictionary.Add(counter++,
                            this.controlPointsDictionary[j* this.NumberOfControlPointsZeta+
                            i* this.NumberOfControlPointsZeta* this.NumberOfControlPointsHeta]);
                    }
                }
                FacesDictionary.Add(2, faceBottom);
                #endregion
                #region FaceUp
                Face faceUp = new Face();
                faceUp.Degrees[0] = this.DegreeKsi;
                faceUp.Degrees[1] = this.DegreeHeta;
                faceUp.KnotValueVectors.Add(0, this.KnotValueVectorKsi);
                faceUp.KnotValueVectors.Add(1, this.KnotValueVectorHeta);
                faceUp.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    for (int j = 0; j < this.NumberOfControlPointsHeta; j++)
                    {
                        faceUp.ControlPointsDictionary.Add(counter++,
                            this.controlPointsDictionary[this.NumberOfControlPointsZeta-1+j * this.NumberOfControlPointsZeta +
                            i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta]);
                    }
                }
                FacesDictionary.Add(3, faceUp);
                #endregion
                #region FaceFront
                Face faceFront = new Face();
                faceFront.Degrees[0] = this.DegreeKsi;
                faceFront.Degrees[1] = this.DegreeZeta;
                faceFront.KnotValueVectors.Add(0, this.KnotValueVectorKsi);
                faceFront.KnotValueVectors.Add(1, this.KnotValueVectorZeta);
                faceFront.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
                    {
                        faceFront.ControlPointsDictionary.Add(counter++,
                            this.controlPointsDictionary[j+i* this.NumberOfControlPointsHeta* this.NumberOfControlPointsZeta]);
                    }
                }
                FacesDictionary.Add(4, faceFront);
                #endregion
                #region FaceBack
                Face faceBack = new Face();
                faceBack.Degrees[0] = this.DegreeKsi;
                faceBack.Degrees[1] = this.DegreeZeta;
                faceBack.KnotValueVectors.Add(0, this.KnotValueVectorKsi);
                faceBack.KnotValueVectors.Add(1, this.KnotValueVectorZeta);
                faceBack.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
                    {
                        faceBack.ControlPointsDictionary.Add(counter++,
                            this.controlPointsDictionary[j + i * this.NumberOfControlPointsHeta* this.NumberOfControlPointsZeta +
                            this.NumberOfControlPointsZeta*(this.NumberOfControlPointsHeta-1)]);
                    }
                }
                FacesDictionary.Add(5, faceBack);
                #endregion
            }
        }

		internal void CreateNurbsShell()
		{
			BuildEdgesDictionary();
			CreateNURBSShells();
		}

		public void BuildEdgesDictionary()
        {
            if (this.NumberOfDimensions==2)
            {
                #region EdgeRight
                Edge edgeRight = new Edge();
                edgeRight.ID = 0;
                edgeRight.Degree = this.DegreeHeta;
                edgeRight.KnotValueVector = this.KnotValueVectorHeta;
                edgeRight.numberOfControlPoints = this.NumberOfControlPointsHeta;
                edgeRight.Patch = this;
                int counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    edgeRight.ControlPointsDictionary.Add(counter++, this.ControlPointsDictionary[i]);
                }
                EdgesDictionary.Add(0, edgeRight);
                #endregion

                #region EdgeLeft
                Edge edgeLeft = new Edge();
                edgeLeft.ID = 1;
                edgeLeft.Degree = this.DegreeHeta;
                edgeLeft.KnotValueVector = this.KnotValueVectorHeta;
                edgeLeft.numberOfControlPoints = this.NumberOfControlPointsHeta;
                edgeLeft.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    edgeLeft.ControlPointsDictionary.Add(counter++, this.ControlPointsDictionary[i+ this.NumberOfControlPointsHeta*(this.NumberOfControlPointsKsi-1)]);
                }
                EdgesDictionary.Add(1, edgeLeft);
                #endregion

                #region EdgeBottom
                Edge edgeBottom = new Edge();
                edgeBottom.ID = 2;
                edgeBottom.Degree = this.DegreeKsi;
                edgeBottom.KnotValueVector = this.KnotValueVectorKsi;
                edgeBottom.numberOfControlPoints = this.NumberOfControlPointsKsi;
                edgeBottom.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    edgeBottom.ControlPointsDictionary.Add(counter++, this.ControlPointsDictionary[i* this.NumberOfControlPointsHeta]);
                }
                EdgesDictionary.Add(2, edgeBottom);
                #endregion

                #region EdgeUp
                Edge edgeUp = new Edge();
                edgeUp.ID = 3;
                edgeUp.Degree = this.DegreeKsi;
                edgeUp.KnotValueVector = this.KnotValueVectorKsi;
                edgeUp.numberOfControlPoints = this.NumberOfControlPointsKsi;
                edgeUp.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    edgeUp.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i * this.NumberOfControlPointsHeta+ this.NumberOfControlPointsHeta-1]);
                }
                EdgesDictionary.Add(3, edgeUp);
                #endregion
            }
            else
            {
                #region Edge1
                Edge edge1 = new Edge();
                edge1.ID = 0;
                edge1.Degree = this.DegreeZeta;
                edge1.KnotValueVector = this.KnotValueVectorZeta;
                edge1.numberOfControlPoints = this.NumberOfControlPointsZeta;
                edge1.Patch = this;
                int counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
                {
                    edge1.ControlPointsDictionary.Add(counter++, this.ControlPointsDictionary[i]);
                }
                EdgesDictionary.Add(0, edge1);
                #endregion
                #region Edge2
                Edge edge2 = new Edge();
                edge2.ID = 1;
                edge2.Degree = this.DegreeZeta;
                edge2.KnotValueVector = this.KnotValueVectorZeta;
                edge2.numberOfControlPoints = this.NumberOfControlPointsZeta;
                edge2.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
                {
                    edge2.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i+ this.NumberOfControlPointsZeta*(this.NumberOfControlPointsHeta-1)]);
                }
                EdgesDictionary.Add(1, edge2);
                #endregion
                #region Edge3
                Edge edge3 = new Edge();
                edge3.ID = 2;
                edge3.Degree = this.DegreeHeta;
                edge3.KnotValueVector = this.KnotValueVectorHeta;
                edge3.numberOfControlPoints = this.NumberOfControlPointsHeta;
                edge3.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    edge3.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i* this.NumberOfControlPointsZeta]);
                }
                EdgesDictionary.Add(2, edge3);
                #endregion
                #region Edge4
                Edge edge4 = new Edge();
                edge4.ID = 3;
                edge4.Degree = this.DegreeHeta;
                edge4.KnotValueVector = this.KnotValueVectorHeta;
                edge4.numberOfControlPoints = this.NumberOfControlPointsHeta;
                edge4.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    edge4.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i * this.NumberOfControlPointsZeta+ this.NumberOfControlPointsZeta-1]);
                }
                EdgesDictionary.Add(3, edge4);
                #endregion

                #region Edge5
                Edge edge5 = new Edge();
                edge5.ID = 4;
                edge5.Degree = this.DegreeZeta;
                edge5.KnotValueVector = this.KnotValueVectorZeta;
                edge5.numberOfControlPoints = this.NumberOfControlPointsZeta;
                int offset = this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta *
                    (this.NumberOfControlPointsKsi - 1);
                edge5.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
                {
                    edge5.ControlPointsDictionary.Add(counter++, this.ControlPointsDictionary[i+offset]);
                }
                EdgesDictionary.Add(4, edge5);
                #endregion
                #region Edge6
                Edge edge6 = new Edge();
                edge6.ID = 5;
                edge6.Degree = this.DegreeZeta;
                edge6.KnotValueVector = this.KnotValueVectorZeta;
                edge6.numberOfControlPoints = this.NumberOfControlPointsZeta;
                edge6.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
                {
                    edge6.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i + this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)+offset]);
                }
                EdgesDictionary.Add(5, edge6);
                #endregion
                #region Edge7
                Edge edge7 = new Edge();
                edge7.ID = 6;
                edge7.Degree = this.DegreeHeta;
                edge7.KnotValueVector = this.KnotValueVectorHeta;
                edge7.numberOfControlPoints = this.NumberOfControlPointsHeta;
                edge7.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    edge7.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i * this.NumberOfControlPointsZeta+offset]);
                }
                EdgesDictionary.Add(6, edge7);
                #endregion
                #region Edge8
                Edge edge8 = new Edge();
                edge8.ID = 7;
                edge8.Degree = this.DegreeHeta;
                edge8.KnotValueVector = this.KnotValueVectorHeta;
                edge8.numberOfControlPoints = this.NumberOfControlPointsHeta;
                edge8.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
                {
                    edge8.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i * this.NumberOfControlPointsZeta + this.NumberOfControlPointsZeta - 1+offset]);
                }
                EdgesDictionary.Add(7, edge8);
                #endregion

                #region Edge9
                Edge edge9 = new Edge();
                edge9.ID = 8;
                edge9.Degree = this.DegreeKsi;
                edge9.KnotValueVector = this.KnotValueVectorKsi;
                edge9.numberOfControlPoints = this.NumberOfControlPointsKsi;
                edge9.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    edge9.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i* this.NumberOfControlPointsZeta* this.NumberOfControlPointsHeta]);
                }
                EdgesDictionary.Add(8, edge9);
                #endregion
                #region Edge10
                Edge edge10 = new Edge();
                edge10.ID = 9;
                edge10.Degree = this.DegreeKsi;
                edge10.KnotValueVector = this.KnotValueVectorKsi;
                edge10.numberOfControlPoints = this.NumberOfControlPointsKsi;
                edge10.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    edge10.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta
                        + this.NumberOfControlPointsZeta-1]);
                }
                EdgesDictionary.Add(9, edge10);
                #endregion
                #region Edge11
                Edge edge11 = new Edge();
                edge11.ID = 10;
                edge11.Degree = this.DegreeKsi;
                edge11.KnotValueVector = this.KnotValueVectorKsi;
                edge11.numberOfControlPoints = this.NumberOfControlPointsKsi;
                edge11.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    edge11.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta +
                        this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)]);
                }
                EdgesDictionary.Add(10, edge11);
                #endregion
                #region Edge12
                Edge edge12 = new Edge();
                edge12.ID = 11;
                edge12.Degree = this.DegreeKsi;
                edge12.KnotValueVector = this.KnotValueVectorKsi;
                edge12.numberOfControlPoints = this.NumberOfControlPointsKsi;
                edge12.Patch = this;
                counter = 0;
                for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
                {
                    edge12.ControlPointsDictionary.Add(counter++,
                        this.ControlPointsDictionary[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta +
                        this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)+ this.NumberOfControlPointsZeta-1]);
                }
                EdgesDictionary.Add(11, edge12);
                #endregion
            }

        }

        public void BuildControlPointsDictionary()
        {
            List<int> controlPointIDs = new List<int>();
            Dictionary<int, ControlPoint> controlPoints = new Dictionary<int, ControlPoint>();
            foreach (Element element in elementsDictionary.Values)
                foreach (ControlPoint controlPoint in element.ControlPoints)
                {
                    controlPointIDs.Add(controlPoint.ID);
                    if (!controlPoints.ContainsKey(controlPoint.ID))
                        controlPoints.Add(controlPoint.ID, controlPoint);
                }
            controlPointIDs = new List<int>(controlPointIDs.Distinct<int>());
            controlPointIDs.Sort();
            foreach (int controlPointID in controlPointIDs)
                if (!controlPointsDictionary.ContainsKey(controlPointID))
                    controlPointsDictionary.Add(controlPointID, controlPoints[controlPointID]);
        }

        #endregion

        public void ResetMaterialsModifiedProperty()
        {
            this.MaterialsModified = false;
            foreach (Element element in elementsDictionary.Values) element.ElementType.ResetMaterialModified();
        }

        internal void CreatePatchData()
        {
            if (this.NumberOfDimensions == 2)
            {
                CreatePatchData2D();
            }else
            {
                CreatePatchData3D();
            }
        }

        private void CreatePatchData2D()
        {
            CreateNURBSElements2D();
            BuildEdgesDictionary();
        }

        private void CreateNURBSElements2D()
        {
            #region Knots
            Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

            List<Knot> knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = 0.0 });
                    id++;
                }
            }
            #endregion

            #region Elements
            Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
            Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
            if (numberOfElementsKsi * numberOfElementsHeta == 0)
            {
                throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
            }

            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    IList<Knot> knotsOfElement = new List<Knot>();
                    knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j]);
                    knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j + 1]);
                    knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j]);
                    knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j + 1]);

                    int multiplicityElementKsi = 0;
                    if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
                    {
                        multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.DegreeKsi;
                    }

                    int multiplicityElementHeta = 0;
                    if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
                    {
                        multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
                    }

                    int nurbsSupportKsi = this.DegreeKsi + 1;
                    int nurbsSupportHeta = this.DegreeHeta + 1;

                    IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                    for (int k = 0; k < nurbsSupportKsi; k++)
                    {
                        for (int l = 0; l < nurbsSupportHeta; l++)
                        {
                            int controlPointID = (i + multiplicityElementKsi) * this.NumberOfControlPointsHeta +
                                (j + multiplicityElementHeta) + k * this.NumberOfControlPointsHeta + l;
                            elementControlPoints.Add(this.controlPointsDictionary[controlPointID]);
                        }
                    }
                    int elementID = i * numberOfElementsHeta + j;
                    Element element = new NURBSElement2D()
                    {
                        ID = elementID,
                        Patch = this,
                        ElementType = new NURBSElement2D()
                    };
                    element.AddKnots(knotsOfElement);
                    element.AddControlPoints(elementControlPoints);
                    this.elementsDictionary.Add(elementID, element);
                    //this.PatchesDictionary[1].ElementsDictionary.Add(element.ID, element);
                }
            }
            #endregion
        }

		private void CreateNURBSShells()
		{
			#region Knots
			Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
			Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

			List<Knot> knots = new List<Knot>();

			int id = 0;
			for (int i = 0; i < singleKnotValuesKsi.Length; i++)
			{
				for (int j = 0; j < singleKnotValuesHeta.Length; j++)
				{
					knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = 0.0 });
					id++;
				}
			}
			#endregion

			#region Elements
			Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
			Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];

			int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
			int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
			if (numberOfElementsKsi * numberOfElementsHeta == 0)
			{
				throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
			}

			for (int i = 0; i < numberOfElementsKsi; i++)
			{
				for (int j = 0; j < numberOfElementsHeta; j++)
				{
					IList<Knot> knotsOfElement = new List<Knot>();
					knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j]);
					knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j + 1]);
					knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j]);
					knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j + 1]);

					int multiplicityElementKsi = 0;
					if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
					{
						multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.DegreeKsi;
					}

					int multiplicityElementHeta = 0;
					if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
					{
						multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
					}

					int nurbsSupportKsi = this.DegreeKsi + 1;
					int nurbsSupportHeta = this.DegreeHeta + 1;

					IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

					for (int k = 0; k < nurbsSupportKsi; k++)
					{
						for (int l = 0; l < nurbsSupportHeta; l++)
						{
							int controlPointID = (i + multiplicityElementKsi) * this.NumberOfControlPointsHeta +
								(j + multiplicityElementHeta) + k * this.NumberOfControlPointsHeta + l;
							elementControlPoints.Add(this.controlPointsDictionary[controlPointID]);
						}
					}
					int elementID = i * numberOfElementsHeta + j;
					Element element = new NURBSKirchhoffLoveShellElement()
					{
						ID = elementID,
						Patch = this,
						ElementType = new NURBSKirchhoffLoveShellElement()
					};
					element.AddKnots(knotsOfElement);
					element.AddControlPoints(elementControlPoints);
					this.elementsDictionary.Add(elementID, element);
					//this.PatchesDictionary[1].ElementsDictionary.Add(element.ID, element);
				}
			}
			#endregion
		}

		private void CreatePatchData3D()
        {
            CreateNURBSElements3D();
            BuildEdgesDictionary();
            BuildFacesDictionary();
        }

        private void CreateNURBSElements3D()
        {
            #region Knots
            Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];

            List<Knot> knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    for (int k = 0; k < singleKnotValuesZeta.Length; k++)
                    {
                        knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = singleKnotValuesZeta[k] });
                        id++;
                    }

                }
            }
            #endregion

            #region Elements
            Vector singlesKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
            Vector singlesKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];
            Vector singlesKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singlesKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singlesKnotValuesHeta.Length - 1;
            int numberOfElementsZeta = singlesKnotValuesZeta.Length - 1;

            if (numberOfElementsKsi * numberOfElementsHeta * numberOfElementsZeta == 0)
            {
                throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
            }

            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    for (int k = 0; k < numberOfElementsZeta; k++)
                    {
                        IList<Knot> knotsOfElement = new List<Knot>();
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k + 1]);
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k + 1]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k + 1]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k + 1]);

                        int multiplicityElementKsi = 0;
                        if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
                        {
                            multiplicityElementKsi = (int)multiplicityKsi[i + 1] - DegreeKsi;
                        }

                        int multiplicityElementHeta = 0;
                        if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
                        {
                            multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
                        }

                        int multiplicityElementZeta = 0;
                        if (multiplicityZeta[k + 1] - this.DegreeZeta > 0)
                        {
                            multiplicityElementZeta = (int)multiplicityZeta[k + 1] - this.DegreeZeta;
                        }

                        int nurbsSupportKsi = this.DegreeKsi + 1;
                        int nurbsSupportHeta = this.DegreeHeta + 1;
                        int nurbsSupportZeta = this.DegreeZeta + 1;

                        IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                        for (int l = 0; l < nurbsSupportKsi; l++)
                        {
                            for (int m = 0; m < nurbsSupportHeta; m++)
                            {
                                for (int n = 0; n < nurbsSupportZeta; n++)
                                {
                                    int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta * NumberOfControlPointsZeta +
                                        (j + multiplicityElementHeta) * NumberOfControlPointsZeta + (k + multiplicityElementZeta) +
                                        l * NumberOfControlPointsHeta * NumberOfControlPointsZeta + m * NumberOfControlPointsZeta + n;

                                    elementControlPoints.Add(this.controlPointsDictionary[controlPointID]);
                                }
                            }
                        }

                        int elementID = i * numberOfElementsHeta * numberOfElementsZeta + j * numberOfElementsZeta + k;
                        Element element = new NURBSElement3D()
                        {
	                        ID = elementID,
	                        Patch = this,
	                        ElementType = new NURBSElement3D()
                        };
                        element.AddKnots(knotsOfElement);
                        element.AddControlPoints(elementControlPoints);
                        this.elementsDictionary.Add(elementID, element);
                        //this.PatchesDictionary[0].ElementsDictionary.Add(element.ID, element);
                    }
                }
            }
            #endregion
        }

        private Vector FindParametricCoordinates(int degree, IVector knotValueVector)
        {
            if (degree <= 0)
            {
                throw new NotSupportedException("Negative Degree.");
            }
            else if (knotValueVector == null)
            {
                throw new ArgumentNullException("Knot Value Vector is null.");
            }

            int numberOfControlPoints = knotValueVector.Length - degree - 1;

            Vector parametricCoordinates = new Vector(numberOfControlPoints);

            for (int i = 0; i < numberOfControlPoints; i++)
            {
                int leftID = (int)Math.Floor(i + (degree + 1) / 2.0);
                int rightID = (int)Math.Ceiling(i + (degree + 1) / 2.0);

                parametricCoordinates[i] = (knotValueVector[leftID] + knotValueVector[rightID]) / 2.0;
            }

            return parametricCoordinates;

        }

        public double[] GetLocalVectorFromGlobal(Element element, IVector globalVector)
        {
            int localDOFs = 0;
            foreach (IList<DOFType> dofs in element.ElementType.DOFEnumerator.GetDOFTypes(element)) localDOFs += dofs.Count;
            var localVector = new double[localDOFs];

            int pos = 0;
            for (int i=0; i<element.ElementType.DOFEnumerator.GetDOFTypes(element).Count; i++)
            {
                ControlPoint controlPoint = element.ControlPoints[i];
                foreach (DOFType dofType in element.ElementType.DOFEnumerator.GetDOFTypes(element)[i])
                {
                    int dof = ControlPointDOFsDictionary[controlPoint.ID][dofType];
                    if (dof != -1) localVector[pos] = globalVector[dof];
                    pos++;
                }
            }
            return localVector;
        }
    }
}
