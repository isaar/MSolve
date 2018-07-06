using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Materials;

//TODO: Materials should be passed in the constructor, not each method
//TODO: not sure about the quadratures for mass. At least as many Gauss points as there are nodes are necessary 
//      (for positive definite mass matrices), but perhaps I need more for exact integration of the polynomials
//TODO: find an triangular quadrature with 6 points for Tri6 consistent mass. Symmetric triangular quadratures are the most 
//      efficient, but that doesn't concern us for mass matrices. The pdf with symmetric triangular quadratures has some other 
//      ones, as do the course notes: https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf.
namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Creates isoparametric continuum elements with uniform thickness. Abstracts the interpolations, integrations,
    /// extrapolations and any other strategies that differentiate the elements (e.g. Quad4 from Tri6). It is also very 
    /// convenient when the thickness and material properties are the same throughout the whole domain or a region.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ContinuumElement2DFactory
    {
        private static readonly IReadOnlyDictionary<CellType2D, IGaussPointExtrapolation2D> extrapolations;
        private static readonly IReadOnlyDictionary<CellType2D, IQuadrature2D> integrationsForStiffness;
        private static readonly IReadOnlyDictionary<CellType2D, IQuadrature2D> integrationsForMass;
        private static readonly IReadOnlyDictionary<CellType2D, IIsoparametricInterpolation2D> interpolations;

        private ElasticMaterial2D commonMaterial;
        private DynamicMaterial commonDynamicProperties;
        private double commonThickness;

        static ContinuumElement2DFactory()
        {
            // Mass integrations require as many Gauss points as there are nodes, in order for the consisntent mass matrix to be
            // of full rank (and symmetric positive definite)

            // Collections' declarations
            var interpolations = new Dictionary<CellType2D, IIsoparametricInterpolation2D>();
            var integrationsForStiffness = new Dictionary<CellType2D, IQuadrature2D>();
            var integrationsForMass = new Dictionary<CellType2D, IQuadrature2D>();
            var extrapolations = new Dictionary<CellType2D, IGaussPointExtrapolation2D>();

            // Quad4
            interpolations.Add(CellType2D.Quad4, InterpolationQuad4.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Quad4, GaussLegendre2D.Order2x2);
            integrationsForMass.Add(CellType2D.Quad4, GaussLegendre2D.Order2x2);
            extrapolations.Add(CellType2D.Quad4, ExtrapolationGaussLegendre2x2.UniqueInstance);

            // Quad8
            interpolations.Add(CellType2D.Quad8, InterpolationQuad8.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Quad8, GaussLegendre2D.Order3x3);
            integrationsForMass.Add(CellType2D.Quad8, GaussLegendre2D.Order3x3);
            extrapolations.Add(CellType2D.Quad8, ExtrapolationGaussLegendre3x3.UniqueInstance);

            // Quad9
            interpolations.Add(CellType2D.Quad9, InterpolationQuad9.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Quad9, GaussLegendre2D.Order3x3);
            integrationsForMass.Add(CellType2D.Quad9, GaussLegendre2D.Order3x3);
            extrapolations.Add(CellType2D.Quad9, ExtrapolationGaussLegendre3x3.UniqueInstance);

            // Tri3
            interpolations.Add(CellType2D.Tri3, InterpolationTri3.UniqueInstance);
            integrationsForStiffness.Add(CellType2D.Tri3, TriangleQuadratureSymmetricGaussian.Order1Point1);
            integrationsForMass.Add(CellType2D.Tri3, TriangleQuadratureSymmetricGaussian.Order2Points3);
            extrapolations.Add(CellType2D.Tri3, ExtrapolationGaussTriangular1Point.UniqueInstance);

            // Tri 6
            interpolations.Add(CellType2D.Tri6, InterpolationTri6.UniqueInstance);
            // see https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf, p. 24-13, paragraph "options"
            integrationsForStiffness.Add(CellType2D.Tri6, TriangleQuadratureSymmetricGaussian.Order2Points3);
            integrationsForMass.Add(CellType2D.Tri6, TriangleQuadratureSymmetricGaussian.Order4Points6); 
            extrapolations.Add(CellType2D.Tri6, ExtrapolationGaussTriangular3Points.UniqueInstance);

            // Static field assignments
            ContinuumElement2DFactory.interpolations = interpolations;
            ContinuumElement2DFactory.integrationsForStiffness = integrationsForStiffness;
            ContinuumElement2DFactory.integrationsForMass = integrationsForMass;
            ContinuumElement2DFactory.extrapolations = extrapolations;
        }

        public ContinuumElement2DFactory(double commonThickness, ElasticMaterial2D commonMaterial, 
            DynamicMaterial commonDynamicProperties)
        {
            this.commonThickness = commonThickness;
            this.commonMaterial = commonMaterial;
            this.commonDynamicProperties = commonDynamicProperties;
        }

        public ContinuumElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes)
        {
            return CreateElement(cellType, nodes, commonThickness, commonMaterial, commonDynamicProperties);
        }

        public ContinuumElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes, double thickness,
            ElasticMaterial2D material, DynamicMaterial dynamicProperties)
        {
            var materialsAtGaussPoints = new Dictionary<GaussPoint2D, ElasticMaterial2D>();
            foreach (GaussPoint2D gaussPoint in integrationsForStiffness[cellType].IntegrationPoints)
            {
                materialsAtGaussPoints[gaussPoint] = material.Clone();
            }
            return CreateElement(cellType, nodes, thickness, materialsAtGaussPoints, dynamicProperties);
        }

        public ContinuumElement2D CreateElement(CellType2D cellType, IReadOnlyList<Node2D> nodes, double thickness,
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
        {
            //TODO: check if nodes - interpolation and Gauss points - materials match
            return new ContinuumElement2D(thickness, nodes, interpolations[cellType],
                integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType],
                materialsAtGaussPoints, dynamicProperties);
        }
    }
}
