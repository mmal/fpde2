u = {
   "u",
   {1,2},
   { "dirichlet", "abc" }
}

s = {
   type = "solver",
   name = "abc",
   nx = 1,
   mesh = {
      type = "m1",
      ghosted = true
   },
   functions = {u,{"v",3},"g"}
}
