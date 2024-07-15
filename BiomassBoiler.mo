package BiomassBoiler
  package Functions
    function ArrheniusEquation
      input Real A;
      input Real E;
      input Modelica.Units.SI.Temperature T;
      output Real k;
    algorithm
      k := A * exp(-E / (T * 8.314));
    end ArrheniusEquation;

    function 'CO/CO2'
      input Modelica.Units.SI.Temperature T;
      output Real omiga;
    algorithm
      omiga := 2 * (1 + 4.3 * exp(-3390 / T)) / (2 + 4.3 * exp(-3390 / T));
    end 'CO/CO2';
  end Functions;

  package Units
    type MassFraction = Real(final quantity= "MassFraction", final unit="kg/kg", displayUnit="kg/kg", nominal=1, min=0);
  end Units;

  package Basics
    package Interfaces
      connector Fuel_inlet "Fuel inlet connector"
        import Modelica.Units.SI.MassFlowRate;
        import Modelica.Units.SI.Temperature;
        import Modelica.Units.SI.SpecificHeatCapacity;
        import Modelica.Units.SI.SpecificEnergy;
        import Modelica.Units.SI.Density;
        import BiomassBoiler.Units.MassFraction;

        MassFlowRate m_flow "Fuel mass flow rate";
        Temperature T "Fuel temperature";
        SpecificEnergy LHV "Lower heating value";
        SpecificHeatCapacity cp;
        Density rho "Fuel density";
        Modelica.Units.SI.Density rho_bulk "堆积密度";
        MassFraction elements[5] "元素成分{C,H,O,N,S}";
        MassFraction components[4] "组分质量分数，水分、挥发分、固定碳、灰分";

        annotation (Icon(graphics={Rectangle(
                extent={{-100,-100},{100,100}},
                lineColor={0,0,0},
                fillPattern=FillPattern.VerticalCylinder,
                fillColor={0,0,0}), Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={0,128,255})}));
      end Fuel_inlet;

      connector Fuel_outlet
        extends Fuel_inlet;
        annotation (Icon(graphics={Rectangle(
                extent={{-100,-100},{100,100}},
                lineColor={0,0,0},
                fillPattern=FillPattern.VerticalCylinder,
                fillColor={0,0,0}), Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={255,255,255})}));
      end Fuel_outlet;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Interfaces;
  end Basics;

  package Test
    model Model5
      Components.BedCombustion1 bedCombustion1_1 annotation (Placement(
            transformation(origin={-120.000000,-10.000000}, extent={{-10.000000,
                -10.000000},{10.000000,10.000000}})));
      Components.BedCombustion1 bedCombustion1_2 annotation (Placement(
            transformation(origin={-55.00000000000001,-10}, extent={{-10,-10},{
                10,10}})));
      Components.BedCombustion1 bedCombustion1_3 annotation (Placement(
            transformation(origin={10.000000,-10.000000}, extent={{-10.000000,-10.000000},
                {10.000000,10.000000}})));
      Modelica.Blocks.Sources.Constant const[4](k={0.149,0.6797,0.1423,0.029}*2.2075)
        annotation (Placement(transformation(origin = {-200.000000, -10.000000}, extent = {{-10.000000, -10.000000}, {10.000000, 10.000000}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=383.15)
        annotation (Placement(transformation(origin={-119.99999999999999,86},
    extent={{-10,-10},{10,10}},
    rotation=-90)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation(Gr=4.44)
        annotation (Placement(transformation(origin={-119.99999999999999,38},
    extent={{-10,-10},{10,10}},
    rotation=90)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T=405.15)
        annotation (Placement(transformation(origin={-55.000000000000014,86},
    extent={{-10,-10},{10,10}},
    rotation=-90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation1(Gr=4.44)
        annotation (Placement(transformation(origin={-55.000000000000014,38},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature2(T=427.15)
        annotation (Placement(transformation(origin={9,80},
    extent={{-10,-10},{10,10}},
    rotation=-90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation2(Gr=4.44)
        annotation (Placement(transformation(origin={7,42},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
    equation
      connect(bedCombustion1_1.m_sj_out, bedCombustion1_2.m_sj_in)
      annotation(Line(origin={-88,-10},
      points={{-21,0},{21,0}},
      color={0,0,127}));
      connect(const.y, bedCombustion1_1.m_sj_in)
      annotation(Line(origin={-160,-10},
      points={{-29,0},{28,0}},
      color={0,0,127}));
      connect(fixedTemperature.port, bodyRadiation.port_b)
      annotation(Line(origin={-120,62},
      points={{1.4210854715202004e-14,14},{1.4210854715202004e-14,-14}},
      color={191,0,0}));
      connect(fixedTemperature1.port, bodyRadiation1.port_b)
      annotation(Line(origin={-55.00000000000003,62},
    points={{1.4210854715202004e-14,14},{1.4210854715202004e-14,-14}},
    color={191,0,0}));
      connect(bedCombustion1_2.m_sj_out, bedCombustion1_3.m_sj_in)
      annotation(Line(origin={-23,-10},
      points={{-21.000000000000007,0},{21,0}},
      color={0,0,127}));
      connect(bodyRadiation.port_a, bedCombustion1_1.port_a)
      annotation(Line(origin={-120,4},
      points={{1.4210854715202004e-14,24},{1.4210854715202004e-14,-24},{0,-24}},
      color={191,0,0}));
      connect(bodyRadiation1.port_a, bedCombustion1_2.port_a)
      annotation(Line(origin={-55,4},
      points={{-1.4210854715202004e-14,24},{-1.4210854715202004e-14,-24},{-7.105427357601002e-15,-24}},
      color={191,0,0}));
      connect(bodyRadiation2.port_a, bedCombustion1_3.port_a) annotation (Line(
            points={{7,32},{6,32},{6,4},{30,4},{30,-26},{10,-26},{10,-20}},
            color={191,0,0}));
      connect(fixedTemperature2.port, bodyRadiation2.port_b)
        annotation (Line(points={{9,70},{7,70},{7,52}}, color={191,0,0}));
      annotation(Diagram(coordinateSystem(extent={{-100,-100},{100,100}},
    grid={2,2})));
    end Model5;

    model Model4
      import BiomassBoiler.Components.BedCombustion;
      BedCombustion bedCombustion(n=4,
        tao=24,
        T(start=293.15))
        annotation (Placement(transformation(origin = {-80.000000, 10.000000}, extent = {{-10.000000, -10.000000}, {10.000000, 10.000000}})));
      BedCombustion bedCombustion1(n=4,
        tao=24,
        T(start=293.15))
        annotation (Placement(transformation(origin={-36,10},
    extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant const(k=2.2075)
        annotation (Placement(transformation(origin={-184,28},
    extent={{-10,-10},{10,10}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T(
            displayUnit="degC") = 373.15)
        annotation (Placement(transformation(origin={-80,106},
    extent={{-10,-10},{10,10}},
    rotation=-90)));
      BiomassBoiler.Components.BoundaryConditions.FuelSource fuelSource(
        variable_m_flow=true,
        variable_T=false,
        variable_components=false) annotation (Placement(transformation(origin=
                {-122,10}, extent={{-10,-10},{10,10}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature1(T(
            displayUnit="degC") = 473.15)
        annotation (Placement(transformation(origin={-36,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation3(Gr=4.44)
        annotation (Placement(transformation(origin={-36,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation1(Gr=4.44)
        annotation (Placement(transformation(origin={-80,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation(Gr=4.44)
        annotation (Placement(transformation(origin={8,58},
    extent={{-10,-10},{10,10}},
    rotation=90)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature2(T(
            displayUnit="degC") = 573.15)
        annotation (Placement(transformation(origin={8,106},
    extent={{10,-10},{-10,10}},
    rotation=90)));

      BedCombustion bedCombustion2
        annotation (Placement(transformation(origin={8,10},
    extent={{-10,-10},{10,10}})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation2(Gr=4.44)
        annotation (Placement(transformation(origin={52,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature3(T(
            displayUnit="degC") = 673.15)
        annotation (Placement(transformation(origin={52,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Components.BedCombustion bedCombustion3
        annotation (Placement(transformation(origin={52,10.000000000000004},
    extent={{-10,-10},{10,10}})),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation4(Gr=4.44)
        annotation (Placement(transformation(origin={96,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature4(T(
            displayUnit="degC") = 773.15)
        annotation (Placement(transformation(origin={96,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Components.BedCombustion bedCombustion4
        annotation (Placement(transformation(origin={96,10},
    extent={{-10,-10},{10,10}})),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation5(Gr=4.44)
        annotation (Placement(transformation(origin={140,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature5(T(
            displayUnit="degC") = 873.15)
        annotation (Placement(transformation(origin={140,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Components.BedCombustion bedCombustion5
        annotation (Placement(transformation(origin={140,10},
    extent={{-10,-10},{10,10}})),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation6(Gr=4.44)
        annotation (Placement(transformation(origin={184,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature6(T(
            displayUnit="degC") = 973.15)
        annotation (Placement(transformation(origin={184,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Components.BedCombustion bedCombustion6
        annotation (Placement(transformation(origin={184,10},
    extent={{-10,-10},{10,10}})),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation7(Gr=4.44)
        annotation (Placement(transformation(origin={228,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature7(T(
            displayUnit="degC") = 1073.15)
        annotation (Placement(transformation(origin={228,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Components.BedCombustion bedCombustion7
        annotation (Placement(transformation(origin={228,10},
    extent={{-10,-10},{10,10}})),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation8(Gr=4.44)
        annotation (Placement(transformation(origin={272,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature8(T(
            displayUnit="degC") = 1173.15)
        annotation (Placement(transformation(origin={272,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Components.BedCombustion bedCombustion8
        annotation (Placement(transformation(origin={272,10},
    extent={{-10,-10},{10,10}})),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation9(Gr=4.44)
        annotation (Placement(transformation(origin={316,58},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature9(T(
            displayUnit="degC") = 1273.15)
        annotation (Placement(transformation(origin={316,106},
    extent={{10,-10},{-10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Components.BedCombustion bedCombustion9
        annotation (Placement(transformation(origin={316,10},
    extent={{-10,-10},{10,10}})),__MWORKS(BlockSystem(StateMachine)));
    equation






















































































      connect(fuelSource.fuel_outlet, bedCombustion.fuel_in)
      annotation(Line(origin={-110,10},
    points={{-2,0},{19,0}},
    color={0,0,0}));
      connect(bodyRadiation3.port_b, fixedTemperature1.port)
      annotation(Line(origin={39,-20},
    points={{-75,88},{-75,116}},
    color={191,0,0}));
      connect(bodyRadiation3.port_a, bedCombustion1.port_b)
      annotation(Line(origin={-36,34},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(fixedTemperature.port, bodyRadiation1.port_b)
      annotation(Line(origin={-80,78},
    points={{0,18},{0,-10}},
    color={191,0,0}));
      connect(bodyRadiation1.port_a, bedCombustion.port_b)
      annotation(Line(origin={-80,34},
      points={{0,14},{0,-14}},
      color={191,0,0}));
      connect(const.y, fuelSource.m_flow)
      annotation(Line(origin={-152,22},
      points={{-21,6},{20,6},{20,-6}},
      color={0,0,127}));
      connect(bedCombustion.fuel_out, bedCombustion1.fuel_in)
      annotation(Line(origin={-73,10},
    points={{4,0},{26,0}},
    color={0,0,0}));
      connect(bedCombustion1.fuel_out, bedCombustion2.fuel_in)
      annotation(Line(origin={1,10},
    points={{-26,0},{-4,0}},
    color={0,0,0}));
      connect(bodyRadiation.port_a, bedCombustion2.port_b)
      annotation(Line(origin={17,39},
    points={{-9,9},{-9,-19}},
    color={191,0,0}));
      connect(fixedTemperature2.port, bodyRadiation.port_b)
      annotation(Line(origin={8,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bodyRadiation2.port_a, bedCombustion3.port_b)
      annotation(Line(origin={61,39},
    points={{-9,9},{-9,-18.999999999999996}},
    color={191,0,0}));
      connect(fixedTemperature3.port, bodyRadiation2.port_b)
      annotation(Line(origin={52,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bodyRadiation4.port_a, bedCombustion4.port_b)
      annotation(Line(origin={105,39},
    points={{-9,9},{-9,-19}},
    color={191,0,0}));
      connect(fixedTemperature4.port, bodyRadiation4.port_b)
      annotation(Line(origin={96,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bodyRadiation5.port_a, bedCombustion5.port_b)
      annotation(Line(origin={149,39},
    points={{-9,9},{-9,-19}},
    color={191,0,0}));
      connect(fixedTemperature5.port, bodyRadiation5.port_b)
      annotation(Line(origin={140,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bodyRadiation6.port_a, bedCombustion6.port_b)
      annotation(Line(origin={193,39},
    points={{-9,9},{-9,-19}},
    color={191,0,0}));
      connect(fixedTemperature6.port, bodyRadiation6.port_b)
      annotation(Line(origin={184,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bodyRadiation7.port_a, bedCombustion7.port_b)
      annotation(Line(origin={237,39},
    points={{-9,9},{-9,-19}},
    color={191,0,0}));
      connect(fixedTemperature7.port, bodyRadiation7.port_b)
      annotation(Line(origin={228,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bodyRadiation8.port_a, bedCombustion8.port_b)
      annotation(Line(origin={281,39},
    points={{-9,9},{-9,-19}},
    color={191,0,0}));
      connect(fixedTemperature8.port, bodyRadiation8.port_b)
      annotation(Line(origin={272,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bodyRadiation9.port_a, bedCombustion9.port_b)
      annotation(Line(origin={325,39},
    points={{-9,9},{-9,-19}},
    color={191,0,0}));
      connect(fixedTemperature9.port, bodyRadiation9.port_b)
      annotation(Line(origin={316,82},
    points={{0,14},{0,-14}},
    color={191,0,0}));
      connect(bedCombustion2.fuel_out, bedCombustion3.fuel_in)
      annotation(Line(origin={30,10},
      points={{-11,0},{11,0},{11,3.552713678800501e-15}},
      color={0,0,0}));
      connect(bedCombustion3.fuel_out, bedCombustion4.fuel_in)
      annotation(Line(origin={74,10},
      points={{-11,3.552713678800501e-15},{11,3.552713678800501e-15},{11,0}},
      color={0,0,0}));
      connect(bedCombustion4.fuel_out, bedCombustion5.fuel_in)
      annotation(Line(origin={118,10},
      points={{-11,0},{11,0}},
      color={0,0,0}));
      connect(bedCombustion5.fuel_out, bedCombustion6.fuel_in)
      annotation(Line(origin={162,10},
      points={{-11,0},{11,0}},
      color={0,0,0}));
      connect(bedCombustion6.fuel_out, bedCombustion7.fuel_in)
      annotation(Line(origin={206,10},
      points={{-11,0},{11,0}},
      color={0,0,0}));
      connect(bedCombustion7.fuel_out, bedCombustion8.fuel_in)
      annotation(Line(origin={250,10},
      points={{-11,0},{11,0}},
      color={0,0,0}));
      connect(bedCombustion8.fuel_out, bedCombustion9.fuel_in)
      annotation(Line(origin={294,10},
      points={{-11,0},{11,0}},
      color={0,0,0}));
      annotation(Diagram(coordinateSystem(extent={{-100,-100},{100,100}},
    grid={2,2})), experiment(
          StopTime=1000,
          __Dymola_NumberOfIntervals=1000,
          __Dymola_Algorithm="Dassl"));
    end Model4;

    model DynamicFuelBed

      Components.BedCombustion bedCombustion
        annotation (Placement(transformation(extent={{-74,-28},{-54,-8}})));
      Components.BoundaryConditions.FuelSource               fuelSource(
        variable_m_flow=true,
        variable_T=false,
        variable_components=false) annotation (Placement(transformation(origin={-112,-18},
                           extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant const(k=2.2075)
        annotation (Placement(transformation(origin={-162,-16},
    extent={{-10,-10},{10,10}})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation1(Gr=4.44)
        annotation (Placement(transformation(origin={-64,20},
    extent={{-10,-10},{10,10}},
    rotation=90)),__MWORKS(BlockSystem(StateMachine)));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T(
            displayUnit="degC") = 373.15)
        annotation (Placement(transformation(origin={-64,66},
    extent={{-10,-10},{10,10}},
    rotation=-90)));
    equation
      connect(fuelSource.fuel_outlet, bedCombustion.fuel_in)
        annotation (Line(points={{-102,-18},{-75,-18}}, color={0,0,0}));
      connect(const.y, fuelSource.m_flow) annotation (Line(points={{-151,-16},{
              -137,-16},{-137,-12},{-122,-12}}, color={0,0,127}));
      connect(fixedTemperature.port, bodyRadiation1.port_b)
        annotation (Line(points={{-64,56},{-64,30}}, color={191,0,0}));
      connect(bodyRadiation1.port_a, bedCombustion.port_b)
        annotation (Line(points={{-64,10},{-64,-8}}, color={191,0,0}));
    end DynamicFuelBed;

    model A
      Real x;
      Real y;

      function f1
        input Real a;
        output Real b;
      algorithm
        b := cos(a);
      end f1;

      function der_f1 = der(f1, a);

    equation
      x = Modelica.Constants.pi / 2;
      y = der_f1(x);
      annotation (experiment(Interval=1, __Dymola_Algorithm="Dassl"));
    end A;

    model Unnamed
      Components.BedCombustion bedCombustion(m_0=10)
        annotation (Placement(transformation(extent={{-44,24},{-24,44}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T=423.15)
        annotation (Placement(transformation(extent={{-142,28},{-122,48}})));
      Components.BoundaryConditions.FuelSource               fuelSource(
        variable_m_flow=true,
        variable_T=false,
        variable_components=false) annotation (Placement(transformation(origin={-106,-4},
                           extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant const(k=2.2075)
        annotation (Placement(transformation(origin={-168,14},
    extent={{-10,-10},{10,10}})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation(Gr=
            4.44)
        annotation (Placement(transformation(extent={{-96,28},{-76,48}})));
    equation
    //   bedCombustion.C = heatCapacitorWithTwoPort.C;
    //   bedCombustion.C = heatCapacitorWithTwoPort1.C;
    //   thermalConductor.G = bedCombustion.G;
    //   thermalConductor2.G = bedCombustion.G;
      connect(const.y, fuelSource.m_flow) annotation (Line(points={{-157,14},{-120,14},
              {-120,2},{-116,2}}, color={0,0,127}));
      connect(fuelSource.fuel_outlet, bedCombustion.fuel_in) annotation (Line(
            points={{-96,-4},{-50,-4},{-50,34},{-45,34}}, color={0,0,0}));
      connect(fixedTemperature.port, bodyRadiation.port_a)
        annotation (Line(points={{-122,38},{-96,38}}, color={191,0,0}));
      connect(bodyRadiation.port_b, bedCombustion.port_b) annotation (Line(
            points={{-76,38},{-50,38},{-50,44},{-34,44}}, color={191,0,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Unnamed;
  end Test;

  package Components
    model BedCombustion1
      //import BiomassBoiler.Components.MassFlow.Interfaces;
      // import Modelica.SIunits.MassFlowRate;
      // import Modelica.SIunits.Mass;
      import Modelica.Units.SI.*;
      import Modelica.Blocks.Interfaces;
      import BiomassBoiler.Functions.ArrheniusEquation;

      extends BiomassBoiler.Components.HeatCapacities.HeatCapacitorWithNoEnergyEquation;
      extends Modelica.Blocks.Icons.Block;

      parameter Real tao = 24 "停留时间";
      parameter Real cp = 2695 "燃料比热容";

      parameter Integer n = 4 "Number of inputs (= number of outputs)";
      Mass m(start = 1, fixed = true) "体积内物质总质量";
      MassFlowRate m_in;
      MassFlowRate m_out;
      Mass m_sj[n] "区域各种组分质量，水分、挥发分、固定碳、灰分";
      Interfaces.RealInput m_sj_in[n](unit="kg/s") annotation (Placement(
            transformation(extent={{-140,-20},{-100,20}})));
      Interfaces.RealOutput m_sj_out[n](unit="kg/s") annotation (Placement(
            transformation(extent={{100,-10},{120,10}})));

      Real k_evp "燃料水分析出速率常数(s^-1)";
      Real R_evp "燃料水分蒸发速率";
      Real k_vol "燃料挥发分析出速率常数(s^-1)";
      Real R_vol "燃料挥发分析出速率";
      //Real k_cahr "燃料焦炭燃烧速率常数(s^-1)";
      //Real R_char "燃料焦炭燃烧速率";

    initial equation
      //m = 1;
    equation

      // 蒸发速率、挥发分析出速率、焦炭燃烧速率
      k_evp = ArrheniusEquation(5.13 * 10^10, 88000, T);
      R_evp = k_evp * m_sj[1];
      k_vol = ArrheniusEquation(7 * 10^7, 129766, T);
      R_vol = k_vol * m_sj[2];
      //k_cahr =

      // 组分平衡
      m_sj_out = m_sj / tao;
      der(m_sj[1]) = m_sj_in[1] - m_sj_out[1] - R_evp;
      der(m_sj[2]) = m_sj_in[2] - m_sj_out[2] - R_vol;
      der(m_sj[3]) = m_sj_in[3] - m_sj_out[3];
      der(m_sj[4]) = m_sj_in[4] - m_sj_out[4];

      // 总质量、进出量
      m = sum(m_sj);
      m_in = sum(m_sj_in);
      m_out = sum(m_sj_out);

      C = cp * m "热容计算J/K";
    end BedCombustion1;

    model BedCombustion
      import Modelica.Units.SI.*;
      import Modelica.Blocks.Interfaces;
      import BiomassBoiler.Functions.ArrheniusEquation;
      import BiomassBoiler.Basics.Interfaces.Fuel_inlet;

      extends
        BiomassBoiler.Components.HeatCapacities.HeatCapacitorWithNoEnergyEquation;
      extends Modelica.Blocks.Icons.Block;

      parameter Length width = 3.7 "炉排宽度";
      Length length "炉排单元长度";
      //Real V_fuel;
      parameter Velocity v = 15 "炉排速度m/h";
      parameter Real tao = 24 "停留时间";
      parameter Real cp = 2695 "燃料比热容";
      parameter Mass m_0 = 1 "初始燃料质量";
      parameter Integer n_units = 2 "床层离散数量";
      parameter Real G = 3.4 "热导率";
    //   Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=3.4);
    //   Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=2000);
    //   HeatCapacities.HeatCapacitor heatCapacitor[n_units];
    //   ThermalConductor thermalConductor[n_units];


      BiomassBoiler.Basics.Interfaces.Fuel_inlet fuel_in annotation (Placement(transformation(origin={-110,0},
    extent={{-10,-10},{10,10}})));
      BiomassBoiler.Basics.Interfaces.Fuel_outlet fuel_out annotation (Placement(
      transformation(extent={{100,-10},{120,10}})));
      parameter Integer n = 4 "Number of input、output";
      Mass m_sj[n](start = {0.149, 0.6797, 0.1423, 0.029} * m_0) "区域各种组分质量，水分、挥发分、固定碳、灰分";
      Mass m "体积内物质总质量";
      Height bedHeight "床层高度";

      HeatFlowRate Q_fuel_h "进入燃料吸热";
      HeatFlowRate Q_evp "水蒸发吸热";
      Real k_evp "燃料水分析出速率常数(s^-1)";
      Real R_evp "燃料水分蒸发速率";
      Real k_vol "燃料挥发分析出速率常数(s^-1)";
      Real R_vol "燃料挥发分析出速率";
      //Real k_cahr "燃料焦炭燃烧速率常数(s^-1)";
      //Real R_char "燃料焦炭燃烧速率";


    protected
      constant HeatCapacity LHW = 2260000 "水的汽化潜热";
      MassFlowRate m_sj_in[n];
      MassFlowRate m_sj_out[n];

    initial equation

    equation

      /* Fuel inlet */
      m_sj_in = fuel_in.m_flow * fuel_in.components;

      /* Fuel outlet */
      fuel_out.m_flow = sum(m_sj_out);
      fuel_out.components = m_sj / m;
      fuel_out.T = T;
    //   if m <= 0 then
    //     fuel_out.components = zeros(4);
    //   else
    //     fuel_out.components = m_sj / m;
    //   end if;
    //   fuel_out.T = T;

      // 蒸发速率、挥发分析出速率、焦炭燃烧速率
      //k_evp = ArrheniusEquation(5.13 * 10^10, 88000, T);
      k_evp = ArrheniusEquation(4.5 * 10^3, 45000, T);
      R_evp = k_evp * m_sj[1];
      k_vol = ArrheniusEquation(7 * 10^7, 129766, T);
      R_vol = k_vol * m_sj[2];
      //k_cahr =

      // 组分平衡
      m_sj_out = m_sj / tao;
      der(m_sj[1]) = m_sj_in[1] - m_sj_out[1] - R_evp;
      der(m_sj[2]) = m_sj_in[2] - m_sj_out[2] - R_vol;
      der(m_sj[3]) = m_sj_in[3] - m_sj_out[3];
      der(m_sj[4]) = m_sj_in[4] - m_sj_out[4];

      // 能量平衡
      C = cp * m / (n_units + 1) "热容计算J/K";
      //Q_fuel_h = cp * fuel_in.m_flow * (fuel_in.T - T) / tao;
      Q_fuel_h = 0;
      Q_evp = -R_evp * LHW;
      C * der(T) = Q_fuel_h + port_a.Q_flow + port_b.Q_flow + Q_evp;
      //C * der(T) = Q_fuel_h + port_a.Q_flow + Q_evp;
    //   heatCapacitor.C = C;
      // 床层高度间热传导
      //G = 0.2 * width *  length / bedHeight;

    //   for i in 1:n_units loop
    //     thermalConductor[i].G = G;
    //     heatCapacitor[i].C = C;
    //   end for;

    //   connect(port_a, thermalConductor.port_a);
    //   connect(heatCapacitor.port, thermalConductor.port_b);
    //   connect(port_a, thermalConductor[1].port_a);
    //   connect(heatCapacitor[1].port, thermalConductor[1].port_b);
    //   for i in 2:n_units loop
    //     connect(heatCapacitor[i-1].port, thermalConductor[i].port_a);
    //     connect(heatCapacitor[i].port, thermalConductor[i].port_b);
    //   end for;

      // 总质量、床层高度、单元长度、燃料体积
      m = sum(m_sj);
      length = tao * v / 3600;
      bedHeight = m / (width * length) / 500;
      //V_fuel = length * bedHeight * width;


    end BedCombustion;

    package HeatCapacities
      model HeatCapacitorWithNoEnergyEquation
        "Lumped thermal element storing heat"

        Modelica.Units.SI.HeatCapacity C "Heat capacity of element (= cp*m)";
        Modelica.Units.SI.Temperature T(start=293.15, displayUnit="degC")
          "Temperature of element";
        Modelica.Units.SI.TemperatureSlope der_T(start=0)
          "Time derivative of temperature (= der(T))";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation (
            Placement(transformation(
              origin={0,-100},
              extent={{-10,-10},{10,10}},
              rotation=90)));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation (
            Placement(transformation(origin={0,100}, extent={{-10,-10},{10,10}})));
      equation
        T = port_a.T;
        der_T = der(T);
        //C*der(T) = port_a.Q_flow;
        // C*der(T) = port_a.Q_flow + port_b.Q_flow;
        port_b.T = T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={
              Text(
                extent={{-150,110},{150,70}},
                textString="%name",
                lineColor={0,0,255}),
              Polygon(
                points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,13},{
                    -76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},
                    {-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{
                    42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{70,-51},{72,-35},{76,
                    -21},{78,-13},{78,3},{74,15},{66,25},{54,33},{44,41},{36,57},{26,65},
                    {0,67}},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},
                    {-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},
                    {8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,-77},{40,-77},{30,
                    -79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,-73},{-22,-73},{-30,-71},
                    {-40,-65},{-50,-55},{-56,-43},{-58,-35},{-58,-25},{-60,-13},{-60,-5},
                    {-60,7},{-58,17},{-56,19},{-52,27},{-48,35},{-44,45},{-40,57},{-58,
                    35}},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid),
              Text(extent={{-69,7},{71,-24}}, textString="%C")}), Diagram(
              coordinateSystem(
              extent={{-100,-100},{100,100}},
              preserveAspectRatio=true,
              grid={2,2}), graphics={
              Polygon(
                origin={0,-11},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid,
                points={{0,78},{-20,74},{-40,68},{-52,54},{-58,46},{-68,36},{-72,24},{
                    -76,10},{-78,-4},{-76,-20},{-76,-32},{-76,-42},{-70,-54},{-64,-62},
                    {-48,-66},{-30,-72},{-18,-72},{-2,-74},{8,-78},{22,-78},{32,-76},{
                    42,-70},{54,-64},{56,-62},{66,-50},{68,-42},{70,-40},{72,-24},{76,
                    -10},{78,-2},{78,14},{74,26},{66,36},{54,44},{44,52},{36,68},{26,76},
                    {0,78}}),
              Polygon(
                origin={-12,-16},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid,
                points={{-46,51},{-56,41},{-60,29},{-64,15},{-66,1},{-64,-15},{-64,-27},
                    {-64,-37},{-58,-49},{-52,-57},{-36,-61},{-18,-67},{-6,-67},{10,-69},
                    {20,-73},{34,-73},{44,-71},{54,-65},{66,-59},{54,-61},{52,-61},{42,
                    -63},{32,-65},{30,-65},{22,-65},{14,-61},{0,-57},{-10,-57},{-18,-55},
                    {-28,-49},{-38,-39},{-44,-27},{-46,-19},{-46,-9},{-48,3},{-48,11},
                    {-48,23},{-46,33},{-44,35},{-40,43},{-36,51},{-32,61},{-28,73},{-46,
                    51}}),
              Ellipse(
                origin={0,-6.5},
                lineColor={255,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid,
                extent={{-6,5.5},{6,-5.5}}),
              Text(
                origin={30.5,-6},
                extent={{-19.5,19},{19.5,-19}},
                textString="T"),
              Line(
                origin={0,-54},
                points={{0,42},{0,-42}},
                color={255,0,0})}));
      end HeatCapacitorWithNoEnergyEquation;

      model HeatCapacitorWithTwoPort
        Modelica.Units.SI.HeatCapacity C "Heat capacity of element (= cp*m)";
        Modelica.Units.SI.Temperature T(start=293.15, displayUnit="degC")
          "Temperature of element";
        Modelica.Units.SI.TemperatureSlope der_T(start=0)
          "Time derivative of temperature (= der(T))";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a annotation (
            Placement(transformation(
              origin={0,-100},
              extent={{-10,-10},{10,10}},
              rotation=90)));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b annotation (
            Placement(transformation(origin={0,100}, extent={{-10,-10},{10,10}})));
      equation
        T = port_a.T;
        der_T = der(T);
        C*der(T) = port_a.Q_flow + port_b.Q_flow;
        port_b.T = T;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={
              Text(
                extent={{-150,110},{150,70}},
                textString="%name",
                lineColor={0,0,255}),
              Polygon(
                points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,13},{
                    -76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},
                    {-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{
                    42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{70,-51},{72,-35},{76,
                    -21},{78,-13},{78,3},{74,15},{66,25},{54,33},{44,41},{36,57},{26,65},
                    {0,67}},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},
                    {-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},
                    {8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,-77},{40,-77},{30,
                    -79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,-73},{-22,-73},{-30,-71},
                    {-40,-65},{-50,-55},{-56,-43},{-58,-35},{-58,-25},{-60,-13},{-60,-5},
                    {-60,7},{-58,17},{-56,19},{-52,27},{-48,35},{-44,45},{-40,57},{-58,
                    35}},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid),
              Text(extent={{-69,7},{71,-24}}, textString="%C")}), Diagram(
              coordinateSystem(
              extent={{-100,-100},{100,100}},
              preserveAspectRatio=true,
              grid={2,2}), graphics={
              Polygon(
                origin={0,-11},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid,
                points={{0,78},{-20,74},{-40,68},{-52,54},{-58,46},{-68,36},{-72,24},{
                    -76,10},{-78,-4},{-76,-20},{-76,-32},{-76,-42},{-70,-54},{-64,-62},
                    {-48,-66},{-30,-72},{-18,-72},{-2,-74},{8,-78},{22,-78},{32,-76},{
                    42,-70},{54,-64},{56,-62},{66,-50},{68,-42},{70,-40},{72,-24},{76,
                    -10},{78,-2},{78,14},{74,26},{66,36},{54,44},{44,52},{36,68},{26,76},
                    {0,78}}),
              Polygon(
                origin={-12,-16},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid,
                points={{-46,51},{-56,41},{-60,29},{-64,15},{-66,1},{-64,-15},{-64,-27},
                    {-64,-37},{-58,-49},{-52,-57},{-36,-61},{-18,-67},{-6,-67},{10,-69},
                    {20,-73},{34,-73},{44,-71},{54,-65},{66,-59},{54,-61},{52,-61},{42,
                    -63},{32,-65},{30,-65},{22,-65},{14,-61},{0,-57},{-10,-57},{-18,-55},
                    {-28,-49},{-38,-39},{-44,-27},{-46,-19},{-46,-9},{-48,3},{-48,11},
                    {-48,23},{-46,33},{-44,35},{-40,43},{-36,51},{-32,61},{-28,73},{-46,
                    51}}),
              Ellipse(
                origin={0,-6.5},
                lineColor={255,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid,
                extent={{-6,5.5},{6,-5.5}}),
              Text(
                origin={30.5,-6},
                extent={{-19.5,19},{19.5,-19}},
                textString="T"),
              Line(
                origin={0,-54},
                points={{0,42},{0,-42}},
                color={255,0,0})}));
      end HeatCapacitorWithTwoPort;

      model HeatCapacitor "Lumped thermal element storing heat"
        Modelica.Units.SI.HeatCapacity C "Heat capacity of element (= cp*m)";
        Modelica.Units.SI.Temperature T(start=293.15, displayUnit="degC")
          "Temperature of element";
        Modelica.Units.SI.TemperatureSlope der_T(start=0)
          "Time derivative of temperature (= der(T))";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port annotation (
            Placement(transformation(
              origin={0,-100},
              extent={{-10,-10},{10,10}},
              rotation=90)));
      equation
        T = port.T;
        der_T = der(T);
        C*der(T) = port.Q_flow;
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={
              Text(
                extent={{-150,110},{150,70}},
                textString="%name",
                textColor={0,0,255}),
              Polygon(
                points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,13},{
                    -76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},{-64,-73},
                    {-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,-89},{32,-87},{
                    42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{70,-51},{72,-35},{76,
                    -21},{78,-13},{78,3},{74,15},{66,25},{54,33},{44,41},{36,57},{26,65},
                    {0,67}},
                lineColor={160,160,164},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},
                    {-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},
                    {8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,-77},{40,-77},{30,
                    -79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,-73},{-22,-73},{-30,-71},
                    {-40,-65},{-50,-55},{-56,-43},{-58,-35},{-58,-25},{-60,-13},{-60,-5},
                    {-60,7},{-58,17},{-56,19},{-52,27},{-48,35},{-44,45},{-40,57},{-58,
                    35}},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid),
              Text(extent={{-69,7},{71,-24}}, textString="%C")}));
      end HeatCapacitor;
    end HeatCapacities;

    package BoundaryConditions
      model FuelSource
        import Modelica.Units.SI.MassFlowRate;
        import Modelica.Units.SI.Temperature;
        import BiomassBoiler.Units.MassFraction;
        import BiomassBoiler.Basics.Interfaces.*;

        parameter Boolean variable_m_flow = false "True, if mass flow defined by variable input" annotation(Dialog(group="Define Variable Boundaries"));
        parameter Boolean variable_T = false "True, if temperature defined by variable input" annotation(Dialog(group="Define Variable Boundaries"));
        parameter Boolean variable_components = false "True, if composition defined by variable input"    annotation(Dialog(group="Define Variable Boundaries"));

        parameter MassFlowRate m_flow_const = 2.2075 "Constant mass flow rate" annotation (Dialog(group="Constant Boundaries", enable=not variable_m_flow));
        parameter Temperature T_const = 293.15 "Constant specific temperature of source" annotation (Dialog(group="Constant Boundaries", enable=not variable_T));
        parameter MassFraction components_const[4] = {0.149, 0.6797, 0.1423, 0.029} "Constant composition" annotation (Dialog(group="Constant Boundaries", enable=not variable_xi));

      protected
        MassFlowRate m_flow_in;
        Temperature T_in;
        MassFraction components_in[4];

      public
        Fuel_outlet fuel_outlet
          annotation (Placement(transformation(extent={{90,-10},{110,10}})));

        Modelica.Blocks.Interfaces.RealInput m_flow = m_flow_in if (variable_m_flow) "Variable mass flow rate"
          annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
        Modelica.Blocks.Interfaces.RealInput T = T_in if (variable_T) "Variable specific temperature"
          annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
        Modelica.Blocks.Interfaces.RealInput components[4] = components_in
          if (variable_components) "Variable components"
          annotation (Placement(transformation(extent={{-120,-80},{-80,-40}})));

      equation

        if (not variable_m_flow) then
          m_flow_in = m_flow_const;
        end if;
        if (not variable_T) then
          T_in = T_const;
        end if;
        if (not variable_components) then
          components_in = components_const;
        end if;

        fuel_outlet.T = T_in;
        fuel_outlet.m_flow = m_flow_in;
        fuel_outlet.components = components_in;

          annotation (Diagram(graphics={
              Line(points={{40,0},{90,0},{72,10}}),
              Line(points={{90,0},{72,-10}}),
              Rectangle(
                extent={{-40,40},{40,-40}},
                lineColor={0,0,255},
                fillColor={255,255,0},
                fillPattern=FillPattern.CrossDiag),
              Text(extent={{-30,60},{-12,40}}, textString=
                                                   "Q"),
              Text(
                extent={{-64,26},{-40,6}},
                lineColor={0,0,255},
                textString=
                     "P")}),Icon(coordinateSystem(extent={{-100,-100},{100,100}},
      grid={2,2}),graphics={  Rectangle(origin={0,0},
      fillColor={255,255,0},
      fillPattern=FillPattern.Solid,
      extent={{-100,100},{100,-100}})}));
      end FuelSource;
    end BoundaryConditions;

    model BedUnits
      import BiomassBoiler.Components.HeatCapacities.HeatCapacitorWithTwoPort;

      HeatCapacitorWithTwoPort heatCapacitor[10];
    equation

    end BedUnits;

    model ThermalConductor
      "Lumped thermal element transporting heat without storing it"
      extends Modelica.Thermal.HeatTransfer.Interfaces.Element1D;
      Modelica.Units.SI.ThermalConductance G "Constant thermal conductance of material";

    equation
      Q_flow = G*dT;
      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                {100,100}}), graphics={
            Rectangle(
              extent={{-90,70},{90,-70}},
              pattern=LinePattern.None,
              fillColor={192,192,192},
              fillPattern=FillPattern.Backward),
            Line(points={{-90,70},{-90,-70}}, thickness=0.5),
            Line(points={{90,70},{90,-70}}, thickness=0.5),
            Text(
              extent={{-150,120},{150,80}},
              textString="%name",
              textColor={0,0,255}),
            Text(extent={{-150,-80},{150,-110}}, textString="G=%G")}));
    end ThermalConductor;
  end Components;
  annotation (uses(Modelica(version="4.0.0")));
end BiomassBoiler;
