using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Materials;

namespace ISAAR.MSolve.IGA.Readers
{
    public class IsogeometricReader
    {
        enum Attributes
        {
            numberofdimensions,
            numberofpatches,numberofinterfaces,
            patchid,interfaceid,
            degreeksi,degreeheta,degreezeta,
            numberofcpksi,numberofcpheta,numberofcpzeta,
            knotvaluevectorksi,knotvaluevectorheta, knotvaluevectorzeta,
            patchcpid,
            cpcoord,thickness,material,end
        }

        public ModelCreator ModelCreator { get; private set; }
        public string Filename { get; private set; }

        public IsogeometricReader(ModelCreator modelCreator, string filename)
        {
            ModelCreator = modelCreator;
            Filename = filename;
        }

        private int patchID=-1;
        private int numberOfValues=0;
        int[] localControlPointIDs;


        public void CreateModelFromFile()
        {
            char[] delimeters = { ' ', '=', '\t' };
            Attributes? name = null;

            String[] text = System.IO.File.ReadAllLines(Filename);

            for (int i=0; i<text.Length; i++)
            {
                String[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
                if (line.Length == 0)
                {
                    continue;
                }
                try
                {
                    name = (Attributes)Enum.Parse(typeof(Attributes), line[0].ToLower());
                }catch (Exception exception)
                {
                    throw new KeyNotFoundException("Variable name " + line[0] + " is not found.");
                }

                switch (name)
                {
                    case Attributes.numberofdimensions:
                        ModelCreator.NumberOfDimensions = Int32.Parse(line[1]);
                        break;
                    case Attributes.thickness:
                        ModelCreator.Thickness = Double.Parse(line[1], CultureInfo.InvariantCulture);
                        break;
                    case Attributes.numberofpatches:
                        ModelCreator.NumberOfPatches = Int32.Parse(line[1]);
                        break;
                    case Attributes.numberofinterfaces:
                        ModelCreator.Model.NumberOfInterfaces = Int32.Parse(line[1]);
                        break;

                    case Attributes.material:
                        if (ModelCreator.NumberOfDimensions == 2)
                        {
                            ModelCreator.Material = new ElasticMaterial2D((line[4] == "plstress") ? StressState2D.PlaneStress : StressState2D.PlaneStrain) { YoungModulus = Double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = Double.Parse(line[3], CultureInfo.InvariantCulture)};
                        }
                        else
                        {
                            ModelCreator.Material = new ElasticMaterial3D { YoungModulus = Double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = Double.Parse(line[3], CultureInfo.InvariantCulture) };
                        }
                        break;
                    case Attributes.patchid:
                        patchID = Int32.Parse(line[1]);
                        break;

                    case Attributes.degreeksi:
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("Degree Ksi of a patch must be defined after the patchID");
                        ModelCreator.DegreeKsiDictionary.Add(patchID, Int32.Parse(line[1]));
                        break;
                    case Attributes.degreeheta:
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("Degree Heta of a patch must be defined after the patchID");
                        ModelCreator.DegreeHetaDictionary.Add(patchID, Int32.Parse(line[1]));
                        break;
                    case Attributes.degreezeta:
                        if (ModelCreator.NumberOfDimensions == 2)
                            throw new ArgumentOutOfRangeException("You must not define degree Zeta in case of 2D");
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("Degree Zeta of a patch must be defined after the patchID");
                        ModelCreator.DegreeZetaDictionary.Add(patchID, Int32.Parse(line[1]));
                        break;
                    case Attributes.numberofcpksi:
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("Number of Control Points Ksi of a patch must be defined after the patchID");
                        ModelCreator.NumberOfControlPointsKsiDictionary.Add(patchID, Int32.Parse(line[1]));
                        break;
                    case Attributes.numberofcpheta:
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("Number of Control Points Heta of a patch must be defined after the patchID");
                        ModelCreator.NumberOfControlPointsHetaDictionary.Add(patchID, Int32.Parse(line[1]));
                        break;
                    case Attributes.numberofcpzeta:
                        if (ModelCreator.NumberOfDimensions == 2)
                            throw new ArgumentOutOfRangeException("You must not define number of Control Points Zeta in case of 2D");
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("Number of Control Points Zeta of a patch must be defined after the patchID");
                        ModelCreator.NumberOfControlPointsZetaDictionary.Add(patchID, Int32.Parse(line[1]));
                        break;
                    case Attributes.knotvaluevectorksi:
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("KnotValue Vector Ksi of a patch must be defined after the patchID");
                        if (ModelCreator.DegreeKsiDictionary[patchID] == 0 || ModelCreator.NumberOfControlPointsKsiDictionary[patchID] == 0)
                        {
                            throw new ArgumentOutOfRangeException("Degree Ksi and number of Control Points per axis Ksi must be defined before Knot Value Vector Ksi.");
                        }
                        numberOfValues = ModelCreator.NumberOfControlPointsKsiDictionary[patchID] + ModelCreator.DegreeKsiDictionary[patchID] + 1;
                        double[] KnotValueVectorKsi = new double[numberOfValues];
                        for (int j = 0; j < numberOfValues; j++)
                        {
                            KnotValueVectorKsi[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
                        }
                        ModelCreator.KnotValueVectorsKsiDictionary.Add(patchID, KnotValueVectorKsi);
                        break;
                    case Attributes.knotvaluevectorheta:
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("KnotValue Vector Heta of a patch must be defined after the patchID");
                        if (ModelCreator.DegreeHetaDictionary[patchID] == 0 || ModelCreator.NumberOfControlPointsHetaDictionary[patchID] == 0)
                        {
                            throw new ArgumentOutOfRangeException("Degree Heta and number of Control Points per axis Heta must be defined before Knot Value Vector Heta.");
                        }
                        numberOfValues = ModelCreator.NumberOfControlPointsHetaDictionary[patchID] + ModelCreator.DegreeHetaDictionary[patchID] + 1;
                        double[] KnotValueVectorHeta = new double[numberOfValues];
                        for (int j = 0; j < numberOfValues; j++)
                        {
                            KnotValueVectorHeta[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
                        }
                        ModelCreator.KnotValueVectorsHetaDictionary.Add(patchID, KnotValueVectorHeta);
                        break;
                    case Attributes.knotvaluevectorzeta:
                        if (ModelCreator.NumberOfDimensions == 2)
                            throw new ArgumentOutOfRangeException("You must not define Knot Value Vector Zeta in case of 2D");
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("KnotValue Vector Zeta of a patch must be defined after the patchID");
                        if (ModelCreator.DegreeZetaDictionary[patchID] == 0 || ModelCreator.NumberOfControlPointsZetaDictionary[patchID] == 0)
                        {
                            throw new ArgumentOutOfRangeException("Degree Zeta and number of Control Points per axis Zeta must be defined before Knot Value Vector Zeta.");
                        }
                        numberOfValues = ModelCreator.NumberOfControlPointsZetaDictionary[patchID] + ModelCreator.DegreeZetaDictionary[patchID] + 1;
                        double[] KnotValueVectorZeta = new double[numberOfValues];
                        for (int j = 0; j < numberOfValues; j++)
                        {
                            KnotValueVectorZeta[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
                        }
                        ModelCreator.KnotValueVectorsZetaDictionary.Add(patchID, KnotValueVectorZeta);
                        break;

                    case Attributes.patchcpid:
                        if (patchID == -1)
                            throw new ArgumentOutOfRangeException("Control Points ID of a patch must be defined after the patchID");
                        int numberOfPatchCP =(ModelCreator.NumberOfDimensions==3)? ModelCreator.NumberOfControlPointsKsiDictionary[patchID] *
                            ModelCreator.NumberOfControlPointsHetaDictionary[patchID] *
                            ModelCreator.NumberOfControlPointsZetaDictionary[patchID]: 
                            ModelCreator.NumberOfControlPointsKsiDictionary[patchID] *
                            ModelCreator.NumberOfControlPointsHetaDictionary[patchID];
                        localControlPointIDs = new int[numberOfPatchCP];
                        for (int j = 0; j < numberOfPatchCP; j++)
                        {
                            localControlPointIDs[j]= Int32.Parse(line[j + 1]);
                        }
                        ModelCreator.ControlPointIDsDictionary.Add(patchID, localControlPointIDs);
                        break;
                    case Attributes.cpcoord:
                        ModelCreator.NumberOfControlPoints = Int32.Parse(line[1]);
                        for (int j = 0; j < ModelCreator.NumberOfControlPoints; j++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int controlPointGlobalID = Int32.Parse(line[0]);
                            double x = Double.Parse(line[1], CultureInfo.InvariantCulture);
                            double y = Double.Parse(line[2], CultureInfo.InvariantCulture);
                            double z = Double.Parse(line[3], CultureInfo.InvariantCulture);
                            double w = Double.Parse(line[4], CultureInfo.InvariantCulture);
                            ControlPoint controlPoint = new ControlPoint()
                            { ID = controlPointGlobalID, X = x, Y = y, Z = z, WeightFactor = w };
                            ModelCreator.ControlPointsDictionary.Add(controlPointGlobalID, controlPoint);
                        }
                        break;
                    case Attributes.end:
                        ModelCreator.CreateModelData();
                        return;
                }
            }            
        }        


    }
}
