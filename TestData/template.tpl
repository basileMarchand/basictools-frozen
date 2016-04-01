****calcul {calcul}

***init_dof_value TP uniform {Ti}

***mesh 
**file {mesh}

***linear_solver {solver}


***python
**script {python_script}

***resolution
**sequence {sequence}
*time {time}
*increment {increment}
*iteration {iteration}
*ratio {ratio}
*algorithm {algorithm}

***table
{table}  

***bc
{bc}


***material
*file this_file 1

****return

***behavior thermal
**conductivity isotropic
k function {conductivity}

**coefficient 
capacity function {coefficient}

***return
