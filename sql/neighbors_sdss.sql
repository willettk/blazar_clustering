CREATE TABLE #UPLOAD(
   up_ra FLOAT,
   up_dec FLOAT,
   up_id int
)
INSERT INTO #UPLOAD
SELECT RA AS UP_RA,DEC AS UP_DEC,search_id AS UP_ID
FROM mydb.all_blazars 
CREATE TABLE #tmp (
			  up_id int,
			   objid bigint
)
INSERT INTO #tmp
EXEC spgetneighbors 10
INSERT INTO mydb.all_blazars_neighbors
select a.*,t.objid as matched_id from #tmp t, mydb.all_blazars a  where t.up_id = a.search_id 
