Class Opm.Production Extends Ens.Production
{

XData ProductionDefinition
{
<Production Name="Opm.Production" TestingEnabled="true" LogGeneralTraceEvents="true">
  <Description></Description>
  <ActorPoolSize>2</ActorPoolSize>
  <Item Name="Python.bp.CompareProcess" Category="" ClassName="Python.bp.CompareProcess" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bp.SdfProcess" Category="" ClassName="Python.bp.SdfProcess" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bp.SmilesProcess" Category="" ClassName="Python.bp.SmilesProcess" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bopka.PkaPredictorOperation" Category="" ClassName="Python.bopka.PkaPredictorOperation" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bordkit.RDKitOperation" Category="" ClassName="Python.bordkit.RDKitOperation" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bosdf.SdfOperation" Category="" ClassName="Python.bosdf.SdfOperation" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bomisc.IUPACOperation" Category="" ClassName="Python.bomisc.IUPACOperation" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bomisc.GenerateImageOperation" Category="" ClassName="Python.bomisc.GenerateImageOperation" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
  <Item Name="Python.bopersist.Persist" Category="" ClassName="Python.bopersist.Persist" PoolSize="1" Enabled="true" Foreground="false" Comment="" LogTraceEvents="false" Schedule="">
  </Item>
</Production>
}

}
