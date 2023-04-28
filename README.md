# MoM_Basics

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://.github.io/MoM_Basics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://deltaeecs.github.io/MoM_Basics.jl/dev/)
[![Build Status](https://github.com/deltaeecs/MoM_Basics.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/deltaeecs/MoM_Basics.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/deltaeecs/MoM_Basics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/deltaeecs/MoM_Basics.jl)


#### 介绍
提供 MoM 包的基础信息，包含几何、网格相关类型和函数的定义、接口，基函数相关信息接口。MoM 本身被拆分为几个独立的包以方便开发时快读编译， 同时避免在无图形化界面使用时调入绘图相关包而导致报错。