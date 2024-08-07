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

    function Cp_char
      input Modelica.Units.SI.Temperature T;
      output Modelica.Units.SI.SpecificHeatCapacity Cp;
    algorithm
      Cp := 1390 + 0.36*T;
    end Cp_char;

    function Cp_ash
      input Modelica.Units.SI.Temperature T;
      output Modelica.Units.SI.SpecificHeatCapacity Cp;
    algorithm
      Cp := 754 - 0.586*(T - 273.15);
    end Cp_ash;

    function Cp_biomass
      input BiomassBoiler.Units.MassFraction components[4];
      input Modelica.Units.SI.Temperature T;
      output Modelica.Units.SI.SpecificHeatCapacity Cp;
    algorithm
      Cp := 4180*components[1] + 2400*components[2] + Cp_char(T)*components[3] +
        Cp_ash(T)*components[4];
    end Cp_biomass;

    function Gas_viscosity
      input Modelica.Units.SI.Temperature T;
      output Real gasViscosity;
    algorithm
      //gasViscosity := 1.98*10^(-5)*(T/300)^(2/3);
      gasViscosity := 4.847*10^(-7)*T^0.64487;
    end Gas_viscosity;

    function Gas_heatConduction
    //   input Modelica.Units.SI.Temperature T;
    //   //output Modelica.Units.SI.Conductivity
    // algorithm
    //
    end Gas_heatConduction;

    function DiffusionCoef_O2
      input Modelica.Units.SI.Temperature T;
      output Modelica.Units.SI.DiffusionCoefficient D;
    algorithm
      D := 0.207*10^(-4)*(T/300)^1.75;
    end DiffusionCoef_O2;

    function R_mix "床层气体混合速率"
      import BiomassBoiler.Units.ConcentrationRate;
      import SIUnits = Modelica.Units.SI;

      input SIUnits.DiffusionCoefficient D_g "挥发分气体扩散系数";
      input SIUnits.Diameter d_p "燃料颗粒直径";
      input SIUnits.Concentration C_i "燃烧气体i浓度";
      input SIUnits.Concentration C_o2;
      input SIUnits.Velocity u_g "气体表观速度";
      input Real epsilon "床层孔隙率";
      input Real omega_i "反应当量系数";
      input Real omega_o2;

      output ConcentrationRate R_mix;
    algorithm
      R_mix := 0.83*(150*D_g*(1 - epsilon)^(2/3)/(d_p^2*epsilon) + 1.75*u_g*(1 -
        epsilon)^(1/3)/(d_p*epsilon))*min(C_i/omega_i, C_o2/omega_o2);
    end R_mix;

    function DiffusionCoef
      import SIUnits = Modelica.Units.SI;

      input SIUnits.Temperature T;
      input SIUnits.Temperature T_ref;
      input SIUnits.Pressure p;
      input SIUnits.Pressure p_ref;
      input SIUnits.DiffusionCoefficient D_ref;
      output SIUnits.DiffusionCoefficient D;
    algorithm
      D := D_ref*(T/T_ref)^1.75*(p_ref/p);
    end DiffusionCoef;

    function Weighted_Average
      input Real x[:] "weight";
      input Real y[:] "value";
      output Real wa;
    algorithm
      wa := sum(x*y);
    end Weighted_Average;
  end Functions;

  package Units
    type MassFraction = Real(final quantity= "MassFraction", final unit="kg/kg", displayUnit="kg/kg", nominal=1, min=0);
    type ConcentrationRate = Real (final quantity="ConcentrationRate", final unit=
            "mol/(m3.s)");
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

      connector FlueGas_inlet "进入烟气"
        import Modelica.Units.SI.MassFlowRate;
        import Modelica.Units.SI.Temperature;
        import Modelica.Units.SI.SpecificHeatCapacity;
        import Modelica.Units.SI.SpecificEnergy;
        import Modelica.Units.SI.Density;
        import BiomassBoiler.Units.MassFraction;

        MassFlowRate m_flow "flueGas mass flow rate";
        Temperature T "flueGas temperature";
        MassFraction composition[6] "烟气组成";

        annotation (   Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}),
                         graphics={
              Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={118,106,98},
                lineThickness=0.5,
                fillColor={118,106,98},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-60,60},{60,-60}},
                lineColor={27,36,42},
                fillColor={27,36,42},
                fillPattern=FillPattern.Solid)}), Diagram(coordinateSystem(
                preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
                                                          graphics));
      end FlueGas_inlet;

      connector FlueGas_outlet
        extends FlueGas_inlet;
        annotation (Icon(graphics={
              Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={118,106,98},
                lineThickness=0.5,
                fillColor={118,106,98},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-80,80},{80,-80}},
                lineColor={118,106,98},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-60,60},{60,-60}},
                lineColor={27,36,42},
                fillColor={27,36,42},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-45,45},{45,-45}},
                lineColor={27,36,42},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}), Diagram(graphics));
      end FlueGas_outlet;

      connector Air_inlet
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Air_inlet;

      connector Air_outlet
        annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
              coordinateSystem(preserveAspectRatio=false)));
      end Air_outlet;
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
      import Modelica.Units.SI;
      import BiomassBoiler.Units.MassFraction;
      import Air = Modelica.Media.Air.ReferenceAir.Air_pT;

      Air.BaseProperties medium(
        T,
        p);

      SI.Temperature T "Temperature";
      MassFraction X[6] "Composition";
      MassFraction X1[4] "Composition";
      SI.Pressure p "Pressure";
      BiomassBoiler.Components.FlueGasObject flueGas(p=p,T=T,X=X);
      BiomassBoiler.Components.AirObject air(p=p,T=T,X=X1);

      Air.Density rho = Air.density(medium.state);

    equation
      T = 273.15;
      X = {0.1383,0.032,0.0688,1 - 0.1383 - 0.032 - 0.0688 - 0.0000000001 - 0.0000000001,
        0.0000000001,0.0000000001};
      X1 = {0.7808, 0.2095, 0.0093, 1-0.7808-0.2095-0.0093};
      p = 1.0133e5;
      der(medium.T) = 0;
      medium.p = p;



    end Unnamed;

    model ExhaustGas
      import Modelica.Media.IdealGases;
      import Modelica.Media.IdealGases.Common;
      import Modelica.Units.SI;
      import BiomassBoiler.Units.MassFraction;
      package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa (
        mediumName="ExhaustGas",
        data={Common.SingleGasesData.O2,Common.SingleGasesData.CO2,Common.SingleGasesData.H2O,
            Common.SingleGasesData.N2,Common.SingleGasesData.Ar,Common.SingleGasesData.SO2},
        fluidConstants={Common.FluidData.O2,Common.FluidData.CO2,Common.FluidData.H2O,
            Common.FluidData.N2,Common.FluidData.Ar,Common.FluidData.SO2},
        substanceNames={"Oxygen","Carbondioxide","Water","Nitrogen","Argon","Sulfurdioxide"},
        reference_X={0.2095, 0.0003, 0.013, 1 - 0.2095 - 0.0003 - 0.0093 - 0.000001 - 0.013, 0.0093, 0.000001});

      SI.Temperature T "Temperature";
      MassFraction X[6] "Composition";
      SI.Pressure p "Pressure";

      Medium.BaseProperties medium(
        T(start=273.15, fixed=true),
        X(start={0.2095, 0.0003, 0.013, 1 - 0.2095 - 0.0003 - 0.0093 - 0.000001 - 0.013, 0.0093, 0.000001}),
        p(start=1.0133e5, fixed=true));

      // 计算混合物的热力学属性
      Medium.SpecificHeatCapacity cp = Medium.specificHeatCapacityCp(medium.state);
      //Modelica.SIunits.SpecificHeatCapacity cv;
      Modelica.Units.SI.Density rho = Medium.density(medium.state);
      Medium.ThermalConductivity thermalConductivity = Medium.thermalConductivity(medium.state);
      Medium.DynamicViscosity mu = Medium.dynamicViscosity(medium.state);
      Medium.PrandtlNumber pr = Medium.prandtlNumber(medium.state);
      Medium.MolarMass mms = Medium.molarMass(medium.state);

    //algorithm
      // 使用MixtureGasNasa组件的函数计算热力学属性
      //cp := Medium.specificHeatCapacityCp(mixture.state);
      //cv := mixture.specificHeatCapacityCv(T, p);
      //rho := mixture.density(T, p);
      //eta := mixture.dynamicViscosity(T, p);
      //lambda := mixture.thermalConductivity(T, p);
    //   package Medium = IdealGases.SingleGases.CO "Medium model";
    //     Medium.BaseProperties medium(
    //     T(start=273.15, fixed=true),
    //     p(start=1.0e5, fixed=true));
    //   Medium.SpecificHeatCapacity cp=Medium.specificHeatCapacityCp(medium.state);
    equation
    //   der(medium.T) = 0;
    //   medium.p = 1.0133e5;
    //   medium.X = {0.1383,0.032,0.0688,1 - 0.1383 - 0.032 - 0.0688 - 0.0000000001 -
    //         0.0000000001,0.0000000001,0.0000000001};
      der(T) = 0;
      p = 1.0133e5;
      X = {0.1383,0.032,0.0688,1 - 0.1383 - 0.032 - 0.0688 - 0.0000000001 -
            0.0000000001,0.0000000001,0.0000000001};
      medium.T = T;
      medium.p = p;
      medium.X = X;
    end ExhaustGas;

    model Unnamed1
      import Modelica.Media.IdealGases;
      import Modelica.Media.IdealGases.Common;
      package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa (
        mediumName="ExhaustGas",
        data={Common.SingleGasesData.CO, Common.SingleGasesData.CO2, Common.SingleGasesData.H2,
         Common.SingleGasesData.CH4, Common.SingleGasesData.C2H6, Common.SingleGasesData.H2O},
        fluidConstants={Common.FluidData.O2, Common.FluidData.CO2, Common.FluidData.H2,
         Common.FluidData.CH4, Common.FluidData.C2H6, Common.FluidData.H2O},
        substanceNames={"CO", "CO2", "H2", "CH4", "C2H6", "H2O"},
        reference_X={0.507123,0.143,0.0045098,0.1537,0.0606672,0.131});

      Medium.BaseProperties medium(
        T(start=600, fixed=true),
        X(start={0.507123,0.143,0.0045098,0.1537,0.0606672,0.131}),
        p(start=1.0133e5, fixed=true));

      // 计算混合物的热力学属性
      Medium.SpecificHeatCapacity cp = Medium.specificHeatCapacityCp(medium.state);
      Modelica.Units.SI.Density rho = Medium.density(medium.state);
      Medium.ThermalConductivity thermalConductivity = Medium.thermalConductivity(medium.state);
      Medium.DynamicViscosity mu = Medium.dynamicViscosity(medium.state);
      Medium.PrandtlNumber pr = Medium.prandtlNumber(medium.state);

    equation
      der(medium.T) = 0;
      medium.p = 1.0133e5;
      medium.X = {0.507123,0.143,0.0045098,0.1537,0.0606672,0.131};


    end Unnamed1;

    model UseExternalMedia
      //   package Medium = ExternalMedia.Media.CoolPropMedium (
      //     mediumName="Water",
      //     substanceNames={"Water"},
      //     ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.ph);
      //   Medium.AbsolutePressure p;
      //   Medium.SpecificEnthalpy h;
      //   Medium.Temperature T;
      //   Medium.Density d;
      // equation
      //   p = 1e5; // 设定压力为1 bar
      //   T = Medium.temperature_ph(p, h);
      //   d = Medium.density_ph(p, h);
      //   // 其他需要的热力学性质也可以类似地获取
      import water = Modelica.Media.Water.WaterIF97_pT;
      import      Modelica.Media.IdealGases.SingleGases.O2;


      water.BaseProperties medium(
        T(start=293.15, fixed=true),
        p(start=1.0133e5, fixed=true));

      // 计算密度
      water.SpecificHeatCapacity cp = water.specificHeatCapacityCp(medium.state);
      Modelica.Units.SI.Density rho = water.density(medium.state);
      water.ThermalConductivity thermalConductivity = water.thermalConductivity(medium.state);
      water.DynamicViscosity mu = water.dynamicViscosity(medium.state);
      water.PrandtlNumber pr = water.prandtlNumber(medium.state);
      water.SpecificEnthalpy h = water.specificEnthalpy(medium.state);
    equation
      der(medium.T) = 0;
      medium.p = 1.0133e5;

    end UseExternalMedia;

    model GasReaction
      import BiomassBoiler.ChemicalReactions.GasPhase.GasSpecies;
      import BiomassBoiler.Functions.DiffusionCoef;
      import BiomassBoiler.Functions.Weighted_Average;
      import BiomassBoiler.Functions.DiffusionCoef_O2;
      import BiomassBoiler.Functions.R_mix;
      import Modelica.Units.SI;
      BiomassBoiler.ChemicalReactions.GasPhase.Solution solution(C(start={10,25,10,10,0,0,0.5}));
      BiomassBoiler.ChemicalReactions.GasPhase.Reaction.'CO + 1/2 * O2 -> CO2' reaction1;
      BiomassBoiler.ChemicalReactions.GasPhase.Reaction.'H2 + 1/2 * O2 -> H2O' reaction2;
      BiomassBoiler.ChemicalReactions.GasPhase.Reaction.'CH4 + 3/2 * O2 -> CO + 2 * H2O' reaction3;
      Real xi[GasSpecies];
      Real D_co;
      Real D_o2;
      Real D_h2;
      Real D_ch4;
      Real D_c2h6;
      Real D_co2;
      Real D_h2o;
      Real D_mix;
    //   Real r_mix[7];
      Real r_mix;
      parameter SI.Pressure p_ref = 1.0133e5;
      parameter SI.Pressure p = 1.0133e5;
      parameter SI.Temperature T_ref = 293.15;
      parameter SI.Temperature T = 1000;
    equation
      connect(reaction1.mixture, solution.mixture);
      connect(reaction2.mixture, solution.mixture);
      connect(reaction3.mixture, solution.mixture);

      D_co = DiffusionCoef(T = T,T_ref = T_ref,p=p,p_ref=p_ref,D_ref= 0.208*10^(-4));
      D_o2 = DiffusionCoef_O2(T);
      D_h2 = DiffusionCoef(T = T,T_ref = T_ref,p=p,p_ref=p_ref,D_ref= 0.756*10^(-4));
      D_ch4 = DiffusionCoef(T = T,T_ref = T_ref,p=p,p_ref=p_ref,D_ref= 0.21*10^(-4));
      D_c2h6 = DiffusionCoef(T = T,T_ref = T_ref,p=p,p_ref=p_ref,D_ref= 0.126*10^(-4));
      D_co2 = DiffusionCoef(T = T,T_ref = T_ref,p=p,p_ref=p_ref,D_ref= 0.16*10^(-4));
      D_h2o = DiffusionCoef(T = T,T_ref = T_ref,p=p,p_ref=p_ref,D_ref= 0.242*10^(-4));
      D_mix = Weighted_Average(xi,{D_co,D_o2,D_h2,D_ch4,D_c2h6,D_co2,D_h2o});
      xi = solution.C / sum(solution.C);
    //   for i in 1:7 loop
    //     r_mix[i] = R_mix(D_mix,0.01,solution.C[i],solution.C[2],0.22,0.4,1,0.5);
    //   end for;
      r_mix =  R_mix(1.70219045*10^(-4),0.02,10,5,0.22,0.4,1,0.5);
      reaction1.R_mix = r_mix;
      reaction2.R_mix = r_mix;
      annotation (experiment(
          StopTime=10,
          __Dymola_NumberOfIntervals=1000,
          Tolerance=1e-12,
          __Dymola_Algorithm="Lsodar"));
    end GasReaction;

    model gasTest
      import Modelica.Units.SI;
      import BiomassBoiler.Units.MassFraction;
      import Modelica.Media.IdealGases;
      import Modelica.Media.IdealGases.Common;


    //   input SI.Temperature T "Temperature";
    //   input MassFraction X[4] "Composition";
    //   input SI.Pressure p "Pressure";

      package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa (
          mediumName="ExhaustGas",
          data={Common.SingleGasesData.N2, Common.SingleGasesData.CO2,
          Common.SingleGasesData.O2, Common.SingleGasesData.H2O},
          fluidConstants={Common.FluidData.N2, Common.FluidData.CO2,
          Common.FluidData.O2,Common.FluidData.H2O},
          substanceNames={"N2", "CO2", "O2", "H2O"});

      Medium.BaseProperties medium(
        T,
        X,
        p);

      // 计算混合物的热力学属性
      Medium.SpecificHeatCapacity cp=Medium.specificHeatCapacityCp(medium.state);
      Medium.Density rho=Medium.density(medium.state);
      Medium.ThermalConductivity thermalConductivity=Medium.thermalConductivity(medium.state);
      Medium.DynamicViscosity mu=Medium.dynamicViscosity(medium.state);
      Medium.PrandtlNumber pr=Medium.prandtlNumber(medium.state);
      Medium.SpecificEnthalpy h = Medium.specificEnthalpy(medium.state);
      //Medium.MolarMass mm = Medium.MolarVolume

    equation
      medium.T = 218+273.15;
      medium.X = {0.6918,0.136,0.0899,0.0823};
      medium.p = 100307.6;
    end gasTest;

    model ReactionTest
      extends Vol;
      import BiomassBoiler.Functions.ArrheniusEquation;
      import Modelica.Units.SI.MolarMass;

    //   type GasSpecies = enumeration(
    //       CO,
    //       O2,
    //       H2,
    //       CH4,
    //       C2H6,
    //       CO2,
    //       H2O);

      Modelica.Units.SI.Concentration C[GasSpecies](min=0) "GasSpecies concentrations";
      Modelica.Units.SI.MolarFlowRate fluegas[GasSpecies](min=0);

    //   constant GasSpecies CO=GasSpecies.CO;
    //   constant GasSpecies O2=GasSpecies.O2;
    //   constant GasSpecies H2=GasSpecies.H2;
    //   constant GasSpecies CH4=GasSpecies.CH4;
    //   constant GasSpecies C2H6=GasSpecies.C2H6;
    //   constant GasSpecies CO2=GasSpecies.CO2;
    //   constant GasSpecies H2O=GasSpecies.H2O;
    protected
      constant MolarMass mm[GasSpecies] = {0.02801,0.032,0.002016,0.016042,0.030068,0.04401,0.018016};
      parameter Modelica.Units.SI.Temperature T=973.15;
      Real k1=ArrheniusEquation(
          3.25*10^7,
          15098*8.314,
          T);
      Real k2=ArrheniusEquation(
          51.8,
          3420*8.314,
          T)*T^1.5;
      Real k3=ArrheniusEquation(
          1.585*10^10,
          24157*8.314,
          T);
      Real k4=ArrheniusEquation(
          2.67*10^8,
          20131*8.314,
          T)*T^0.5;

    public
      Real R_CO;
      Real R_H2;
      Real R_CH4;
      Real R_C2H6;
    initial equation
      C = {10,50,10,10,10,0,0};
    equation

    //   R_CO = k1*C[CO]*C[O2]^0.5*C[H2O]^0.5;
    //   R_H2 = k2*C[H2]^1.5*C[O2];
    //   R_CH4 = k3*C[CH4]^0.7*C[O2]^0.8;
    //   R_C2H6 = k4*C[C2H6]*C[O2];

    //   fluegas = 1.5*y./mm;
      for i in 1:size(fluegas,1) loop
        if i <> 2 then
          fluegas[i] = 1.5*y[i]/mm[i];
        else
          fluegas[i] = (fluegas[CH4]*2 + fluegas[C2H6]*3.5 + fluegas[CO]*0.5 + fluegas[H2]*0.5)*1.2;
        end if;
      end for;
      if noEvent(C[CO] <= 1e-6 or C[O2] <= 1e-6) then
        R_CO = 0;
      else
        R_CO = k1*C[CO]*C[O2]^0.5*C[H2O]^0.5;
      end if;

      if noEvent(C[H2] <= 1e-6 or C[O2] <= 1e-6) then
        R_H2 = 0;
      else
        R_H2 = min(k2*C[H2]^1.5*C[O2], 1279);
      end if;

      if noEvent(C[CH4] <= 1e-6 or C[O2] <= 1e-6) then
        R_CH4 = 0;
      else
        R_CH4 = k3*C[CH4]^0.7*C[O2]^0.8;
      end if;

      if noEvent(C[C2H6] <= 1e-6 or C[O2] <= 1e-6) then
        R_C2H6 = 0;
      else
        R_C2H6 = min(k4*C[C2H6]*C[O2],1279);
      end if;

      der(C[CO]) = -R_CO + R_CH4 + 2*R_C2H6 + fluegas[CO];
      der(C[H2]) = -R_H2 + fluegas[H2];
      der(C[CH4]) = -R_CH4 + fluegas[CH4];
      der(C[O2]) = -R_CO/2 - R_H2/2 - R_CH4*1.5 - R_C2H6*2.5 + fluegas[O2] - C[O2];
      der(C[CO2]) = R_CO + fluegas[CO2] - C[CO2];
      der(C[H2O]) = R_H2 + R_CH4*2 + R_C2H6*3 + fluegas[H2O] - C[H2O];
      der(C[C2H6]) = -R_C2H6 + fluegas[C2H6];

      //   der(C[CO]) = -k1*C[CO]*C[O2]^0.5*C[H2O]^0.5 + k3*C[CH4]^0.7*C[O2]^0.8 + 2*k4*C[C2H6]*C[O2];
      //   der(C[H2]) = -k2*C[H2]^1.5*C[O2];
      //   der(C[CH4]) = -k3*C[CH4]^0.7*C[O2]^0.8;
      //   der(C[O2]) = -k1*C[CO]*C[O2]^0.5*C[H2O]^0.5/2 - k2*C[H2]^1.5*C[O2]/2 - k3*C[CH4]^0.7*C[O2]^0.8*1.5 - 2.5*k4*C[C2H6]*C[O2];
      //   der(C[CO2]) = k1*C[CO]*C[O2]^0.5*C[H2O]^0.5;
      //   der(C[H2O]) = k2*C[H2]^1.5*C[O2] + k3*C[CH4]^0.7*C[O2]^0.8*2 + 3*k4*C[C2H6]*C[O2];
      //   der(C[C2H6]) = -k4*C[C2H6]*C[O2];


      annotation (experiment(
          StopTime=1000,
          Interval=0.2,
          __Dymola_Algorithm="Dassl"));
    end ReactionTest;

    model Vol
      type GasSpecies = enumeration(
          CO,
          O2,
          H2,
          CH4,
          C2H6,
          CO2,
          H2O);
    protected
      constant Real MW_C = 12.01;    // 碳的摩尔质量, kg/kmol
      constant Real MW_H = 1.008;    // 氢的摩尔质量, kg/kmol
      constant Real MW_O = 16.00;    // 氧的摩尔质量, kg/kmol
      constant Real MW_O2 = 32.00;    // O2的摩尔质量, kg/kmol
      constant Real MW_CO = 28.01;   // CO的摩尔质量, kg/kmol
      constant Real MW_CO2 = 44.01;  // CO2的摩尔质量, kg/kmol
      constant Real MW_H2 = 2.016;   // H2的摩尔质量, kg/kmol
      constant Real MW_CH4 = 16.042; // CH4的摩尔质量, kg/kmol
      constant Real MW_C2H6 = 30.068; // C2H6的摩尔质量, kg/kmol
      constant Real MW_H2O = 18.016; // H2O的摩尔质量, kg/kmol

      constant GasSpecies CO=GasSpecies.CO;
      constant GasSpecies O2=GasSpecies.O2;
      constant GasSpecies H2=GasSpecies.H2;
      constant GasSpecies CH4=GasSpecies.CH4;
      constant GasSpecies C2H6=GasSpecies.C2H6;
      constant GasSpecies CO2=GasSpecies.CO2;
      constant GasSpecies H2O=GasSpecies.H2O;
    public
      Modelica.Units.SI.MassFraction y[GasSpecies](min=0);
      Real t;
    equation
      y[O2] = 0;
      y[CO2] = 0.143;
      y[CH4] = 0.1537;
      y[H2O] = 0.131;
    //   y[H2] = 0.02;
      sum(y) = t;
      y[CO]*MW_C/MW_CO + y[CO2]*MW_C/MW_CO2 + y[CH4]*MW_C/MW_CH4 + y[C2H6]*MW_C*2/MW_C2H6 = 0.42;
      y[CO]*MW_O/MW_CO + y[CO2]*MW_O*2/MW_CO2 + y[H2O]*MW_O/MW_H2O = 0.51;
      y[CH4]*MW_H*4/MW_CH4 + y[C2H6]*MW_H*6/MW_C2H6 + y[H2O]*MW_H*2/MW_H2O + y[H2]*MW_H*2/MW_H2 = 0.07;

    end Vol;
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
      Real cp "燃料比热容";
      parameter Mass m_0 = 1 "初始燃料质量";


      BiomassBoiler.Basics.Interfaces.Fuel_inlet fuel_in annotation (Placement(transformation(origin={-110,0},
    extent={{-10,-10},{10,10}})));
      BiomassBoiler.Basics.Interfaces.Fuel_outlet fuel_out annotation (Placement(
      transformation(extent={{100,-10},{120,10}})));
      parameter Integer n = 4 "Number of input、output";
      Mass m_sj[n](start = {0.149, 0.6797, 0.1423, 0.029} * m_0) "区域各种组分质量，水分、挥发分、固定碳、灰分";
      Mass m "体积内物质总质量";
      Height bedHeight "床层高度";
      MassFlowRate m_sj_in[n];
      MassFlowRate m_sj_out[n];

      HeatFlowRate Q_fuel_h "进入燃料吸热";
      HeatFlowRate Q_evp "水蒸发吸热";
      Real k_evp "燃料水分析出速率常数(s^-1)";
      MassFlowRate R_evp "燃料水分蒸发速率";
      Real k_vol "燃料挥发分析出速率常数(s^-1)";
      MassFlowRate R_vol "燃料挥发分析出速率";
      //Real k_cahr "燃料焦炭燃烧速率常数(s^-1)";
      //Real R_char "燃料焦炭燃烧速率";


    protected
      constant HeatCapacity LHW = 1993470 "水的汽化潜热";

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

      // 蒸发速率、挥发分析出速率、焦炭燃烧速率
      //k_evp = ArrheniusEquation(5.13 * 10^10, 88000, T);
      k_evp = ArrheniusEquation(5.6 * 10^8, 88000, T);
      //k_evp = ArrheniusEquation(4.5 * 10^3, 45000, T);
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
      //C = cp * m / (n_units + 1) "热容计算J/K";
      cp = BiomassBoiler.Functions.Cp_biomass(fuel_in.components, fuel_in.T);
      C = cp * m;
      Q_fuel_h = cp * fuel_in.m_flow * (fuel_in.T - T) / tao;
      Q_evp = -R_evp * LHW;
      C * der(T) = Q_fuel_h + port_up.Q_flow + port_down.Q_flow + Q_evp;

      //C * der(T) = Q_fuel_h + port_a.Q_flow + Q_evp;
      // 床层高度间热传导
      //G = 0.2 * width *  length / bedHeight;

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
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_down annotation (
            Placement(transformation(
              origin={0,-100},
              extent={{-10,-10},{10,10}},
              rotation=90)));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_up annotation (
            Placement(transformation(origin={0,100}, extent={{-10,-10},{10,10}})));
      equation
        T = port_down.T;
        der_T = der(T);
        //C*der(T) = port_a.Q_flow;
        // C*der(T) = port_a.Q_flow + port_b.Q_flow;
        port_up.T = T;
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

    model FlueGasObject
      import Modelica.Units.SI;
      import BiomassBoiler.Units.MassFraction;
      import Modelica.Media.IdealGases;
      import Modelica.Media.IdealGases.Common;


      input SI.Temperature T "Temperature";
      input MassFraction X[6] "Composition";
      input SI.Pressure p "Pressure";

      package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa (
          mediumName="ExhaustGas",
          data={Common.SingleGasesData.CO,Common.SingleGasesData.CO2,Common.SingleGasesData.H2,
              Common.SingleGasesData.CH4,Common.SingleGasesData.C2H6,Common.SingleGasesData.H2O},
          fluidConstants={Common.FluidData.O2,Common.FluidData.CO2,Common.FluidData.H2,
              Common.FluidData.CH4,Common.FluidData.C2H6,Common.FluidData.H2O},
          substanceNames={"CO", "CO2", "H2", "CH4", "C2H6", "H2O"});

      Medium.BaseProperties medium(
        T,
        X,
        p);

      // 计算混合物的热力学属性
      Medium.SpecificHeatCapacity cp=Medium.specificHeatCapacityCp(medium.state);
      Medium.Density rho=Medium.density(medium.state);
      Medium.ThermalConductivity thermalConductivity=Medium.thermalConductivity(medium.state);
      Medium.DynamicViscosity mu=Medium.dynamicViscosity(medium.state);
      Medium.PrandtlNumber pr=Medium.prandtlNumber(medium.state);
      //Medium.MolarMass mm = Medium.MolarVolume

    equation
      medium.T = T;
      medium.X = X;
      medium.p = p;
    end FlueGasObject;

    model AirObject
      import Modelica.Units.SI;
      import BiomassBoiler.Units.MassFraction;
      import Modelica.Media.IdealGases;
      import Modelica.Media.IdealGases.Common;


      input SI.Temperature T "Temperature";
      input MassFraction X[4] "Composition";
      input SI.Pressure p "Pressure";

      package Medium = Modelica.Media.IdealGases.Common.MixtureGasNasa (
          mediumName="ExhaustGas",
          data={Common.SingleGasesData.N2,Common.SingleGasesData.O2,Common.SingleGasesData.Ar,
              Common.SingleGasesData.CO2},
          fluidConstants={Common.FluidData.N2,Common.FluidData.O2,Common.FluidData.Ar,
              Common.FluidData.CO2},
          substanceNames={"N2","O2","Ar","CO2"});

      Medium.BaseProperties medium(
        T,
        X,
        p);

      // 计算混合物的热力学属性
      Medium.SpecificHeatCapacity cp=Medium.specificHeatCapacityCp(medium.state);
      Medium.Density rho=Medium.density(medium.state);
      Medium.ThermalConductivity thermalConductivity=Medium.thermalConductivity(
          medium.state);
      Medium.DynamicViscosity mu=Medium.dynamicViscosity(medium.state)
        "动力粘度";
      Medium.PrandtlNumber pr=Medium.prandtlNumber(medium.state) "普朗特数";
      Medium.SpecificEnthalpy h=Medium.specificEnthalpy(medium.state) "比焓";

    equation
      medium.T = T;
      medium.X = X;
      medium.p = p;
    end AirObject;
  end Components;

  package SubSystem
    model Unnamed
      Components.BedCombustion bedCombustion
        annotation (Placement(transformation(extent={{-24,8},{-4,28}})));
      Components.BedCombustion bedCombustion1
        annotation (Placement(transformation(extent={{-20,-66},{0,-46}})));
      Components.BoundaryConditions.FuelSource               fuelSource(
        variable_m_flow=true,
        variable_T=false,
        variable_components=false) annotation (Placement(transformation(origin={-78,-18},
                           extent={{-10,-10},{10,10}})));
      Modelica.Blocks.Sources.Constant const(k=0.22075)
        annotation (Placement(transformation(origin={-140,0},
    extent={{-10,-10},{10,10}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T(
            displayUnit="K") = 1073.15)
        annotation (Placement(transformation(extent={{-96,56},{-76,76}})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation(Gr=0.37)
        annotation (Placement(transformation(extent={{-48,56},{-28,76}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor
        thermalConductor(G=5.55) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-12,-16})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation1(Gr=0.37)
        annotation (Placement(transformation(extent={{10,-24},{30,-4}})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation2(Gr=0.37)
        annotation (Placement(transformation(extent={{-10,-10},{10,10}},
            rotation=270,
            origin={22,-96})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor
        thermalConductor1(G=5.55) annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={-10,-94})));
      Components.BedCombustion bedCombustion2
        annotation (Placement(transformation(extent={{-18,-144},{2,-124}})));
    equation
      connect(const.y, fuelSource.m_flow) annotation (Line(points={{-129,0},{
              -94,0},{-94,-12},{-88,-12}}, color={0,0,127}));
      connect(fuelSource.fuel_outlet, bedCombustion.fuel_in) annotation (Line(
            points={{-68,-18},{-30,-18},{-30,18},{-25,18}}, color={0,0,0}));
      connect(fuelSource.fuel_outlet, bedCombustion1.fuel_in) annotation (Line(
            points={{-68,-18},{-34,-18},{-34,-56},{-21,-56}}, color={0,0,0}));
      connect(bodyRadiation.port_b, bedCombustion.port_b) annotation (Line(
            points={{-28,66},{-14,66},{-14,28}}, color={191,0,0}));
      connect(fixedTemperature.port, bodyRadiation.port_a)
        annotation (Line(points={{-76,66},{-48,66}}, color={191,0,0}));
      connect(bedCombustion.port_a, thermalConductor.port_a)
        annotation (Line(points={{-14,8},{-12,8},{-12,-6}}, color={191,0,0}));
      connect(bedCombustion1.port_b, thermalConductor.port_b) annotation (Line(
            points={{-10,-46},{-12,-46},{-12,-26}}, color={191,0,0}));
      connect(bedCombustion.port_a, bodyRadiation1.port_a) annotation (Line(
            points={{-14,8},{-14,6},{-12,6},{-12,0},{4,0},{4,-14},{10,-14}},
            color={191,0,0}));
      connect(bodyRadiation1.port_b, bedCombustion1.port_b) annotation (Line(
            points={{30,-14},{34,-14},{34,-40},{-10,-40},{-10,-46}}, color={191,
              0,0}));
      connect(bedCombustion1.port_a, thermalConductor1.port_a) annotation (Line(
            points={{-10,-66},{-10,-75},{-10,-75},{-10,-84}}, color={191,0,0}));
      connect(bedCombustion1.port_a, bodyRadiation2.port_a) annotation (Line(
            points={{-10,-66},{-10,-78},{22,-78},{22,-86}}, color={191,0,0}));
      connect(thermalConductor1.port_b, bedCombustion2.port_b) annotation (Line(
            points={{-10,-104},{-8,-104},{-8,-124}}, color={191,0,0}));
      connect(bodyRadiation2.port_b, bedCombustion2.port_b) annotation (Line(
            points={{22,-106},{22,-124},{-8,-124}}, color={191,0,0}));
      connect(fuelSource.fuel_outlet, bedCombustion2.fuel_in) annotation (Line(
            points={{-68,-18},{-34,-18},{-34,-56},{-26,-56},{-26,-120},{-19,
              -120},{-19,-134}}, color={0,0,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Unnamed;

    model BedUnits
      import Modelica.Units.SI.MassFlowRate;
      import Modelica.Units.SI.Mass;
      import Modelica.Units.SI.Length;
      import Modelica.Units.SI.Height;

      parameter Integer n_units = 10;
      //parameter Length l;
      //parameter Length w;
      Height bedHeight;
      MassFlowRate fuel_in_total;
      MassFlowRate fuel_out_total;
      Mass m_total;
      Basics.Interfaces.Fuel_inlet fuel_inlet[n_units]
        annotation (Placement(transformation(extent={{-110,18},{-90,38}})));
      Basics.Interfaces.Fuel_outlet fuel_outlet[n_units]
        annotation (Placement(transformation(extent={{90,18},{110,38}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_down
        annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_up
        annotation (Placement(transformation(extent={{-10,88},{10,108}})));
      Components.BedCombustion bedCombustion[n_units]
        annotation (Placement(transformation(extent={{-10,20},{10,40}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor[n_units-1](G=4.44)
        annotation (Placement(transformation(
            extent={{-10,-10},{10,10}},
            rotation=270,
            origin={0,-16})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation[n_units-1](Gr=0.37)
        annotation (Placement(transformation(extent={{32,-26},{52,-6}})));
    equation
      connect(bedCombustion[1].port_up, port_up);
      connect(bedCombustion[n_units].port_down, port_down);
      for i in 1:n_units-1 loop
        connect(bedCombustion[i].port_down, thermalConductor[i].port_a);
        connect(bedCombustion[i+1].port_up, thermalConductor[i].port_b);
        connect(bedCombustion[i].port_down, bodyRadiation[i].port_a);
        connect(bedCombustion[i+1].port_up, bodyRadiation[i].port_b);
      end for;
      connect(bedCombustion.fuel_in, fuel_inlet);
      connect(bedCombustion.fuel_out, fuel_outlet);
      m_total = sum(bedCombustion.m);
      fuel_in_total = sum(bedCombustion.fuel_in.m_flow);
      fuel_out_total = sum(bedCombustion.fuel_out.m_flow);
      bedHeight = sum(bedCombustion.bedHeight);
    end BedUnits;

    model Unnamed2
      BedUnits unnamed1_1
        annotation (Placement(transformation(extent={{-10,-6},{10,14}})));
      Modelica.Blocks.Sources.Constant const[10](k=0.22075)
        annotation (Placement(transformation(origin={-110,22},
    extent={{-10,-10},{10,10}})));
      Components.BoundaryConditions.FuelSource               fuelSource[10](
        variable_m_flow=true,
        variable_T=false,
        variable_components=false) annotation (Placement(transformation(origin={-48,4},
                           extent={{-10,-10},{10,10}})));
      BedUnits unnamed1_2
        annotation (Placement(transformation(extent={{84,-4},{104,16}})));
      BedUnits unnamed1_3
        annotation (Placement(transformation(extent={{182,-4},{202,16}})));
      BedUnits unnamed1_4
        annotation (Placement(transformation(extent={{284,-2},{304,18}})));
      Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature4(T(
            displayUnit="degC") = 873.15)
        annotation (Placement(transformation(extent={{-60,78},{-40,98}})));
      Modelica.Thermal.HeatTransfer.Components.BodyRadiation bodyRadiation4(Gr=0.37)
        annotation (Placement(transformation(extent={{-12,78},{8,98}})));
      BedUnits unnamed1_5
        annotation (Placement(transformation(extent={{384,0},{404,20}})));
    equation
      connect(const.y, fuelSource.m_flow) annotation (Line(points={{-99,22},{-64,22},
              {-64,10},{-58,10}}, color={0,0,127}));
      connect(fuelSource.fuel_outlet, unnamed1_1.fuel_inlet) annotation (Line(
            points={{-38,4},{-36,4},{-36,6.8},{-10,6.8}}, color={0,0,0}));
      connect(fixedTemperature4.port, bodyRadiation4.port_a)
        annotation (Line(points={{-40,88},{-12,88}}, color={191,0,0}));
      connect(unnamed1_1.fuel_outlet, unnamed1_2.fuel_inlet) annotation (Line(
            points={{10,6.8},{12,6.8},{12,8.8},{84,8.8}}, color={0,0,0}));
      connect(unnamed1_2.fuel_outlet, unnamed1_3.fuel_inlet)
        annotation (Line(points={{104,8.8},{182,8.8}}, color={0,0,0}));
      connect(unnamed1_3.fuel_outlet, unnamed1_4.fuel_inlet) annotation (Line(
            points={{202,8.8},{243,8.8},{243,10.8},{284,10.8}}, color={0,0,0}));
      connect(unnamed1_4.fuel_outlet, unnamed1_5.fuel_inlet) annotation (Line(
            points={{304,10.8},{344,10.8},{344,12.8},{384,12.8}}, color={0,0,0}));
      connect(bodyRadiation4.port_b, unnamed1_1.port_up) annotation (Line(
            points={{8,88},{12,88},{12,18},{0,18},{0,13.8}}, color={191,0,0}));
      connect(bodyRadiation4.port_b, unnamed1_2.port_up) annotation (Line(
            points={{8,88},{12,88},{12,22},{94,22},{94,15.8}}, color={191,0,0}));
      connect(bodyRadiation4.port_b, unnamed1_3.port_up) annotation (Line(
            points={{8,88},{12,88},{12,15.8},{192,15.8}}, color={191,0,0}));
      connect(bodyRadiation4.port_b, unnamed1_4.port_up) annotation (Line(
            points={{8,88},{12,88},{12,17.8},{294,17.8}}, color={191,0,0}));
      connect(bodyRadiation4.port_b, unnamed1_5.port_up) annotation (Line(
            points={{8,88},{12,88},{12,19.8},{394,19.8}}, color={191,0,0}));
      annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
            coordinateSystem(preserveAspectRatio=false)));
    end Unnamed2;
  end SubSystem;

  package ChemicalReactions
    package GasPhase
      type GasSpecies = enumeration(
          CO,
          O2,
          H2,
          CH4,
          C2H6,
          CO2,
          H2O);
      model Solution
        GasPhase.Interfaces.Mixture mixture;
        Modelica.Units.SI.Concentration C[GasSpecies]=mixture.C "CO,O2,H2,CH4,C2H6,CO2,H2O";
      equation
        der(mixture.C) = mixture.R;
      end Solution;

      package Interfaces
        connector Mixture
          Modelica.Units.SI.Concentration C[GasSpecies];
          flow BiomassBoiler.Units.ConcentrationRate R[GasSpecies]
            "concentrationRate";
        end Mixture;

        partial model Reaction
          import BiomassBoiler.Units.ConcentrationRate;
          //parameter Real k "Reaction coefficient";
          Real k "Reaction coefficient";
          Mixture mixture;
        protected
          ConcentrationRate consumed[GasSpecies];
          ConcentrationRate produced[GasSpecies];
          Modelica.Units.SI.Concentration C[GasSpecies]=mixture.C;
        equation
          consumed = -produced;
          mixture.R = consumed;
        end Reaction;
      end Interfaces;

      package Reaction
        model 'CO + 1/2 * O2 -> CO2'
          extends Interfaces.Reaction;
          import BiomassBoiler.Functions.ArrheniusEquation;
          import BiomassBoiler.Units.ConcentrationRate;
        public
          ConcentrationRate R_mix;
        protected
          ConcentrationRate R;
        initial equation
        equation
          k = ArrheniusEquation(
            3.25*10^7,
            15098*8.314,
            1000);
          if noEvent(C[GasSpecies.CO]>0 and C[GasSpecies.O2]>0) then
            R = min(k*C[GasSpecies.CO]*sqrt(C[GasSpecies.O2])*sqrt(C[GasSpecies.H2O]),R_mix);
          else
            R = 0;
          end if;
        //   C[GasSpecies.CO] = noEvent(max(0, C[GasSpecies.CO]));
          consumed[GasSpecies.CO] = R;
          consumed[GasSpecies.O2] = R/2;
          consumed[GasSpecies.H2] = 0;
          consumed[GasSpecies.CH4] = 0;
          consumed[GasSpecies.C2H6] = 0;
          produced[GasSpecies.CO2] = R;
          consumed[GasSpecies.H2O] = 0;

        end 'CO + 1/2 * O2 -> CO2';

        model 'H2 + 1/2 * O2 -> H2O'
          extends Interfaces.Reaction;
          import BiomassBoiler.Functions.ArrheniusEquation;
          import BiomassBoiler.Units.ConcentrationRate;
        public
          ConcentrationRate R_mix;
        protected
          ConcentrationRate R;

        initial equation
        equation
          //   k = ArrheniusEquation(
          //     1.631*10^9,
          //     24157*8.314,
          //     1000);
          k = ArrheniusEquation(
            51.8,
            3420*8.314,
            1000);
          //   if noEvent(C[GasSpecies.H2] > 0 and C[GasSpecies.O2] > 0) then
          //     R = min(k*C[GasSpecies.H2]^1.5*C[GasSpecies.O2]*1000^1.5, R_mix);
          //   else
          //     R = 0;
          //   end if;
          if noEvent(C[GasSpecies.H2] > 0) then
            R = min(k*C[GasSpecies.H2]^1.5*C[GasSpecies.O2]*1000^1.5, R_mix);
          else
            R = 0;
          end if;
          //   R = k*max(C[GasSpecies.H2],0)^1.5*max(C[GasSpecies.O2],0)*1000^1.5;
          //   C[GasSpecies.H2] = smooth(1, max(0,C[GasSpecies.H2]));
          consumed[GasSpecies.CO] = 0;
          //   consumed[GasSpecies.O2] = R;
          //   consumed[GasSpecies.H2] = 2*R;
          consumed[GasSpecies.O2] = R/2;
          consumed[GasSpecies.H2] = R;
          consumed[GasSpecies.CH4] = 0;
          consumed[GasSpecies.C2H6] = 0;
          consumed[GasSpecies.CO2] = 0;
          //   produced[GasSpecies.H2O] = 2*R;
          produced[GasSpecies.H2O] = R;
        end 'H2 + 1/2 * O2 -> H2O';

        model 'CH4 + 3/2 * O2 -> CO + 2 * H2O'
          extends Interfaces.Reaction;
          import BiomassBoiler.Functions.ArrheniusEquation;
          import BiomassBoiler.Units.ConcentrationRate;

        protected
          ConcentrationRate R;
        initial equation
        equation
          k = ArrheniusEquation(
            1.585*10^10,
            24157*8.314,
            1000);
          if noEvent(C[GasSpecies.CH4]>0 and C[GasSpecies.O2]>0) then
            R = k * C[GasSpecies.CH4]^0.7 * C[GasSpecies.O2]^0.8;
          else
            R = 0;
          end if;
          produced[GasSpecies.CO] = R;
          consumed[GasSpecies.O2] = 1.5*R;
          consumed[GasSpecies.H2] = 0;
          consumed[GasSpecies.CH4] = R;
          consumed[GasSpecies.C2H6] = 0;
          consumed[GasSpecies.CO2] = 0;
          produced[GasSpecies.H2O] = 2*R;
        end 'CH4 + 3/2 * O2 -> CO + 2 * H2O';

        model 'C2H6 + 5/2 * O2 -> 2 * CO + 3 * H2O'
          annotation (Icon(coordinateSystem(preserveAspectRatio=false)),
              Diagram(coordinateSystem(preserveAspectRatio=false)));
        end 'C2H6 + 5/2 * O2 -> 2 * CO + 3 * H2O';
      end Reaction;
    end GasPhase;

    package SolidPhase
    end SolidPhase;
  end ChemicalReactions;
  annotation (uses(Modelica(version="4.0.0")));
end BiomassBoiler;
