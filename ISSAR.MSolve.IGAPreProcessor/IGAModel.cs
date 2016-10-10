using ISAAR.MSolve.Matrices;
using ISSAR.MSolve.IGAPreProcessor.Elements;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGAPreProcessor
{
    public class IGAModel
    {
        private int totalDofs;
        private int numberOfDimensions;

        private int numberOfCPKsi;
        private int numberOfCPHeta;
        private int numberOfCPZeta;
        //private int thickness;

        private int degreeKsi;
        private int degreeHeta;
        private int degreeZeta;
        private Vector<double> knotValueVectorKsi;
        private Vector<double> knotValueVectorHeta;
        private Vector<double> knotValueVectorZeta;
        private readonly IList<ControlPoint> controlPoints = new List<ControlPoint>();
        private readonly IList<Knot> knots = new List<Knot>();
        private readonly IList<IGAElement> elements = new List<IGAElement>();
        private readonly IList<IGALoad> loads = new List<IGALoad>();

        #region Properties
        public int TotalDOFs
        {
            get { return this.totalDofs; }
            set { numberOfDimensions=value; }
        }

        public int NumberOfDimensions
        {
            get { return this.numberOfDimensions; }
            set { numberOfDimensions = value; }
        }

        public int DegreeKsi
        {
            get { return this.degreeKsi; }
            set { degreeKsi = value; }
        }

        public int DegreeHeta
        {
            get { return this.degreeHeta; }
            set { degreeHeta = value; }
        }

        public int DegreeZeta
        {
            get { return this.degreeZeta; }
            set { degreeZeta = value; }
        }

        public Vector<double> KnotValueVectorKsi
        {
            get { return this.knotValueVectorKsi; }
            set { knotValueVectorKsi = value; }
        }

        public Vector<double> KnotValueVectorHeta
        {
            get { return this.knotValueVectorHeta; }
            set { knotValueVectorHeta = value; }
        }

        public Vector<double> KnotValueVectorZeta
        {
            get { return this.knotValueVectorZeta; }
            set { knotValueVectorZeta = value; }
        }

        public int NumberOfCPKsi
        {
            get { return numberOfCPKsi; }
            set { numberOfCPKsi = value; }
        }

        public int NumberOfCPHeta
        {
            get { return numberOfCPHeta; }
            set { numberOfCPHeta = value; }
        }

        public int NumberOfCPZeta
        {
            get { return numberOfCPZeta; }
            set { numberOfCPZeta = value; }
        }

        public IList<ControlPoint> ControlPoints
        {
            get { return controlPoints; }
        }

        public IList<Knot> Knots
        {
            get { return knots; }
        }

        public IList<IGAElement> Elements
        {
            get { return elements; }
        }

        public IList<IGALoad> Loads
        {
            get { return loads; }
        }

        #endregion

        #region Data Creation Routines
        public void CreateDataStructures()
        {
            throw new NotImplementedException();

        }
        

        public void CreateModelData(Matrix2D<double> cpCoordinates)
        {
            if (this.numberOfDimensions == 2)
            {
                CreateModelData2D(cpCoordinates);
            }else
            {
                CreateModelData3D(cpCoordinates);
            }
        }

        private void CreateModelData2D(Matrix2D<double> cpCoordinates)
        {
            CreateControlPoints2D(cpCoordinates);
            CreateKnots2D();
            CreateNURBSELements2D();
        }

        private void CreateNURBSELements2D()
        {
            Vector<double> singlesKnotValuesKsi = knotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> multiplicityKsi = knotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
            Vector<double> singlesKnotValuesHeta = knotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> multiplicityHeta = knotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singlesKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singlesKnotValuesHeta.Length - 1;
            if (numberOfElementsKsi * numberOfElementsHeta == 0)
            {
                throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
            }

            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    IList<Knot> knotsOfElement = new List<Knot>();
                    knotsOfElement.Add(knots[i * singlesKnotValuesHeta.Length + j]);
                    knotsOfElement.Add(knots[i * singlesKnotValuesHeta.Length + j + 1]);
                    knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesHeta.Length + j]);
                    knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesHeta.Length + j + 1]);

                    int multiplicityElementKsi = 0;
                    if (multiplicityKsi[i + 1] - this.degreeKsi > 0)
                    {
                        multiplicityElementKsi = (int)multiplicityKsi[i + 1] - degreeKsi;
                    }

                    int multiplicityElementHeta = 0;
                    if (multiplicityHeta[j + 1] - this.degreeHeta > 0)
                    {
                        multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.degreeHeta;
                    }

                    int nurbsSupportKsi = this.degreeKsi + 1;
                    int nurbsSupportHeta = this.degreeHeta + 1;

                    Vector<int> connectivity = new Vector<int>(nurbsSupportKsi * nurbsSupportHeta);
                    int index = 0;

                    for (int k = 0; k < nurbsSupportKsi; k++)
                    {
                        for (int l = 0; l < nurbsSupportHeta; l++)
                        {
                            connectivity[index] = (i + multiplicityElementKsi) * numberOfCPHeta +
                                (j + multiplicityElementHeta) + k * numberOfCPHeta + l;
                            index++;
                        }
                    }
                    int elementID = i * numberOfElementsHeta + j;
                    IGAElement element = new IGAElement() {ElementType= new NURBSElement2D(elementID, knotsOfElement, connectivity)};
                    //NURBSElement2D nurbsElement = new NURBSElement2D(elementID, knotsOfElement, connectivity);
                    this.elements.Add(element);

                }
            }
        }

        private void CreateKnots2D()
        {
            Vector<double> singleKnotValuesKsi = knotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> singleKnotValuesHeta = knotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    Knot knot = new Knot(id, singleKnotValuesKsi[i], singleKnotValuesHeta[j], 0.0);
                    id++;
                    this.knots.Add(knot);
                }
            }
        }

        private void CreateControlPoints2D(Matrix2D<double> cpCoordinates)
        {
            Vector<double> parametricCoordinatesKsi = FindParametricCoordinates(this.degreeKsi, this.knotValueVectorKsi);
            Vector<double> parametricCoordinatesHeta = FindParametricCoordinates(this.degreeHeta, this.knotValueVectorHeta);

            int id = 0;

            for (int i = 0; i < this.numberOfCPKsi; i++)
            {
                for (int j = 0; j < this.numberOfCPHeta; j++)
                {
                    ControlPoint controlPoint = new ControlPoint(id, cpCoordinates[id, 0], cpCoordinates[id, 1], 0.0, parametricCoordinatesKsi[i], parametricCoordinatesHeta[j], 0.0, cpCoordinates[id, 2]);
                    id++;
                    this.controlPoints.Add(controlPoint);
                }
            }
        }

        private void CreateModelData3D(Matrix2D<double> cpCoordinates)
        {
            CreateControlPoints3D(cpCoordinates);
            CreateKnots3D();
            CreateNURBSElements3D();
        }

        private void CreateNURBSElements3D()
        {
            Vector<double> singlesKnotValuesKsi = knotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> multiplicityKsi = knotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
            Vector<double> singlesKnotValuesHeta = knotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> multiplicityHeta = knotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];
            Vector<double> singlesKnotValuesZeta = knotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> multiplicityZeta = knotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[1];

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
                        if (multiplicityKsi[i + 1] - this.degreeKsi > 0)
                        {
                            multiplicityElementKsi = (int)multiplicityKsi[i + 1] - degreeKsi;
                        }

                        int multiplicityElementHeta = 0;
                        if (multiplicityHeta[j + 1] - this.degreeHeta > 0)
                        {
                            multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.degreeHeta;
                        }

                        int multiplicityElementZeta = 0;
                        if (multiplicityZeta[k + 1] - this.degreeZeta > 0)
                        {
                            multiplicityElementZeta = (int)multiplicityZeta[k + 1] - this.degreeZeta;
                        }

                        int nurbsSupportKsi = this.degreeKsi + 1;
                        int nurbsSupportHeta = this.degreeHeta + 1;
                        int nurbsSupportZeta = this.degreeZeta + 1;

                        Vector<int> connectivity = new Vector<int>(nurbsSupportKsi * nurbsSupportHeta * nurbsSupportZeta);
                        int index = 0;

                        for (int l = 0; l < nurbsSupportKsi; l++)
                        {
                            for (int m = 0; m < nurbsSupportHeta; m++)
                            {
                                for (int n = 0; n < nurbsSupportZeta; n++)
                                {
                                    connectivity[index] = (i + multiplicityElementKsi) * numberOfCPHeta * numberOfCPZeta +
                                    (j + multiplicityElementHeta) * numberOfCPZeta + (k + multiplicityElementZeta) +
                                    l * numberOfCPHeta * numberOfCPZeta + m * numberOfCPZeta + n;
                                    index++;
                                }                                
                            }
                        }

                        int elementID = i * numberOfElementsHeta * numberOfElementsZeta + j * numberOfElementsZeta + k;
                        IGAElement element = new IGAElement() { ElementType = new NURBSElement3D(elementID, knotsOfElement, connectivity) };
                        //NURBSElement2D nurbsElement = new NURBSElement3D(elementID, knotsOfElement, connectivity);

                        this.elements.Add(element);
                        
                    }
                }
            }
        }

        private void CreateKnots3D()
        {
            Vector<double> singleKnotValuesKsi = knotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> singleKnotValuesHeta = knotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector<double> singleKnotValuesZeta = knotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    for (int k = 0; k < singleKnotValuesZeta.Length; k++)
                    {
                        Knot knot = new Knot(id, singleKnotValuesKsi[i], singleKnotValuesHeta[j], singleKnotValuesZeta[k]);
                        id++;
                        this.knots.Add(knot);
                    }

                }
            }
        }

        private void CreateControlPoints3D(Matrix2D<double> cpCoordinates)
        {
            Vector<double> parametricCoordinatesKsi = FindParametricCoordinates(this.degreeKsi, this.knotValueVectorKsi);
            Vector<double> parametricCoordinatesHeta = FindParametricCoordinates(this.degreeHeta, this.knotValueVectorHeta);
            Vector<double> parametricCoordinatesZeta = FindParametricCoordinates(this.degreeZeta, this.knotValueVectorZeta);

            int id = 0;

            for (int i = 0; i < this.numberOfCPKsi; i++)
            {
                for (int j = 0; j < this.numberOfCPHeta; j++)
                {
                    for (int k = 0; k < this.numberOfCPZeta; k++)
                    {
                        ControlPoint controlPoint = new ControlPoint(id, cpCoordinates[id, 0], cpCoordinates[id, 1], cpCoordinates[id, 2], parametricCoordinatesKsi[i], parametricCoordinatesHeta[j], parametricCoordinatesZeta[k], cpCoordinates[id, 3]);
                        id++;
                        this.controlPoints.Add(controlPoint);
                    }

                }
            }
        }

        private Vector<double> FindParametricCoordinates(int degree, Vector<double> knotValueVector)
        {
            if (degree <= 0)
            {
                throw new NotSupportedException("Negative Degree.");
            } else if (knotValueVector == null)
            {
                throw new ArgumentNullException("Knot Value Vector is null.");
            }

            int numberOfControlPoints = knotValueVector.Length - degree - 1;

            Vector<double> parametricCoordinates = new Vector<double>(numberOfControlPoints);

            for (int i = 0; i < numberOfControlPoints; i++)
            {
                int leftID = (int)Math.Floor(i + (degree + 1) / 2.0);
                int rightID = (int)Math.Ceiling(i + (degree + 1) / 2.0);

                parametricCoordinates[i] = (knotValueVector[leftID] + knotValueVector[rightID]) / 2.0;
            }

            return parametricCoordinates;

        }

        #endregion

        public void Clear()
        {
            this.controlPoints.Clear();
            this.knots.Clear();
            this.elements.Clear();
        }

    }
}
