//using ISAAR.MSolve.IGA.Entities;
//using ISAAR.MSolve.IGA.Problems.Structural.Constitutive;
//using ISAAR.MSolve.Numerical.LinearAlgebra;
//using System;
//using System.Collections.Generic;
//using System.Globalization;
//using System.Linq;
//using System.Text;
//using System.Threading.Tasks;

//namespace ISAAR.MSolve.IGA.Readers
//{
//    public static class ModelReader
//    {
//        enum Attributes
//        {
//            numberofdimensions,
//            degreeksi, degreeheta, degreezeta,
//            numberofcpksi, numberofcpheta, numberofcpzeta,
//            knotvaluevectorksi, knotvaluevectorheta, knotvaluevectorzeta,
//            cpcoord, end, thickness, material
//        }

//        public static void CreateModelFromFile(Model model, string fileName)
//        {
//            char[] delimeters = { ' ', '=', '\t' };
//            Attributes? name = null;

//            String[] text = System.IO.File.ReadAllLines(fileName);

//            for (int i = 0; i < text.Length; i++)
//            {
//                String[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
//                if (line.Length == 0)
//                {
//                    continue;
//                }
//                try
//                {
//                    name = (Attributes)Enum.Parse(typeof(Attributes), line[0].ToLower());
//                }
//                catch (Exception e)
//                {
//                    throw new KeyNotFoundException("Variable name is not found " + line[0]);
//                }

//                switch (name)
//                {
//                    case Attributes.numberofdimensions:
//                        model.NumberOfDimensions = Int32.Parse(line[1]);
//                        break;
//                    case Attributes.thickness:
//                        model.Thickness = Double.Parse(line[1], CultureInfo.InvariantCulture);
//                        break;
//                    case Attributes.degreeksi:
//                        model.DegreeKsi = Int32.Parse(line[1]);
//                        break;
//                    case Attributes.degreeheta:
//                        model.DegreeHeta = Int32.Parse(line[1]);
//                        break;
//                    case Attributes.degreezeta:
//                        model.DegreeZeta = Int32.Parse(line[1]);
//                        break;
//                    case Attributes.numberofcpksi:
//                        model.NumberOfControlPointsKsi = Int32.Parse(line[1]);
//                        break;
//                    case Attributes.numberofcpheta:
//                        model.NumberOfControlPointsHeta = Int32.Parse(line[1]);
//                        break;
//                    case Attributes.numberofcpzeta:
//                        model.NumberOfControlPointsZeta = Int32.Parse(line[1]);
//                        break;
//                    case Attributes.knotvaluevectorksi:
//                        if (model.DegreeKsi == 0 || model.NumberOfControlPointsKsi == 0)
//                        {
//                            throw new ArgumentOutOfRangeException("Degree Ksi and number of Control Points per axis Ksi must be defined before Knot Value Vector Ksi.");
//                        }
//                        model.KnotValueVectorKsi = new Vector(model.NumberOfControlPointsKsi + model.DegreeKsi + 1);
//                        for (int j = 0; j < model.NumberOfControlPointsKsi + model.DegreeKsi + 1; j++)
//                        {
//                            model.KnotValueVectorKsi[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
//                        }
//                        break;
//                    case Attributes.knotvaluevectorheta:
//                        if (model.DegreeHeta == 0 || model.NumberOfControlPointsHeta == 0)
//                        {
//                            throw new ArgumentOutOfRangeException("Degree Heta and number of Control Points per axis Heta must be defined before Knot Value Vector Heta.");
//                        }
//                        model.KnotValueVectorHeta = new Vector(model.NumberOfControlPointsHeta + model.DegreeHeta + 1);
//                        for (int j = 0; j < model.NumberOfControlPointsHeta + model.DegreeHeta + 1; j++)
//                        {
//                            model.KnotValueVectorHeta[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
//                        }
//                        break;
//                    case Attributes.knotvaluevectorzeta:
//                        if (model.DegreeZeta == 0 || model.NumberOfControlPointsZeta == 0)
//                        {
//                            throw new ArgumentOutOfRangeException("Degree Zeta and number of Control Points per axis Zeta must be defined before Knot Value Vector Zeta.");
//                        }
//                        model.KnotValueVectorZeta = new Vector(model.NumberOfControlPointsZeta + model.DegreeZeta + 1);
//                        for (int j = 0; j < model.NumberOfControlPointsZeta + model.DegreeZeta + 1; j++)
//                        {
//                            model.KnotValueVectorZeta[j] = Double.Parse(line[j + 1], CultureInfo.InvariantCulture);
//                        }
//                        break;
//                    case Attributes.cpcoord:
//                        if (model.NumberOfDimensions == 2)
//                        {
//                            Matrix2D cpCoordinates = new Matrix2D(model.NumberOfControlPointsKsi * model.NumberOfControlPointsHeta, 3);
//                            cpCoordinates[0, 0] = Double.Parse(line[1], CultureInfo.InvariantCulture);
//                            cpCoordinates[0, 1] = Double.Parse(line[2], CultureInfo.InvariantCulture);
//                            cpCoordinates[0, 2] = Double.Parse(line[3], CultureInfo.InvariantCulture);
//                            for (int j = 1; j < model.NumberOfControlPointsKsi * model.NumberOfControlPointsHeta; j++)
//                            {
//                                i++;
//                                line = text[i].Split(delimeters);
//                                cpCoordinates[j, 0] = Double.Parse(line[0], CultureInfo.InvariantCulture);
//                                cpCoordinates[j, 1] = Double.Parse(line[1], CultureInfo.InvariantCulture);
//                                cpCoordinates[j, 2] = Double.Parse(line[2], CultureInfo.InvariantCulture);
//                            }
//                            model.CreateModelData(cpCoordinates);
//                        }
//                        else
//                        {
//                            Matrix2D cpCoordinates = new Matrix2D(model.NumberOfControlPointsKsi * model.NumberOfControlPointsHeta * model.NumberOfControlPointsZeta, 4);
//                            cpCoordinates[0, 0] = Double.Parse(line[1], CultureInfo.InvariantCulture);
//                            cpCoordinates[0, 1] = Double.Parse(line[2], CultureInfo.InvariantCulture);
//                            cpCoordinates[0, 2] = Double.Parse(line[3], CultureInfo.InvariantCulture);
//                            cpCoordinates[0, 3] = Double.Parse(line[4], CultureInfo.InvariantCulture);
//                            for (int j = 1; j < model.NumberOfControlPointsKsi * model.NumberOfControlPointsHeta * model.NumberOfControlPointsZeta; j++)
//                            {
//                                i++;
//                                line = text[i].Split(delimeters);
//                                cpCoordinates[j, 0] = Double.Parse(line[0], CultureInfo.InvariantCulture);
//                                cpCoordinates[j, 1] = Double.Parse(line[1], CultureInfo.InvariantCulture);
//                                cpCoordinates[j, 2] = Double.Parse(line[2], CultureInfo.InvariantCulture);
//                                cpCoordinates[j, 3] = Double.Parse(line[3], CultureInfo.InvariantCulture);
//                            }
//                            model.CreateModelData(cpCoordinates);
//                        }
//                        break;
//                    case Attributes.material:
//                        if (model.NumberOfDimensions == 2)
//                        {
//                            model.Material = new ElasticMaterial2D { YoungModulus = Double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = Double.Parse(line[3], CultureInfo.InvariantCulture), StressState = line[4] };
//                        }else
//                        {
//                            patc.Material = new ElasticMaterial3D { YoungModulus = Double.Parse(line[2], CultureInfo.InvariantCulture), PoissonRatio = Double.Parse(line[3], CultureInfo.InvariantCulture) };
//                        }                        
//                        break;
//                    case Attributes.end:
//                        return;
//                }
//            }
//        }
//    }
//}
