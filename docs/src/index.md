```@meta
DocTestSetup = quote
    using MoM_Basics
end
```

![MoM](./assets/logo.png)
# MoM_Basics

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://deltaeecs.github.io/MoM_Basics.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://deltaeecs.github.io/MoM_Basics.jl/dev/)
[![Build Status](https://github.com/deltaeecs/MoM_Basics.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/deltaeecs/MoM_Basics.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/deltaeecs/MoM_Basics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/deltaeecs/MoM_Basics.jl)

[![CEM_MoMs](https://img.shields.io/badge/CEM_MoMs-github-orange.svg)](https://github.com/deltaeecs/CEM_MoMs.git)
[![CEM_MoMs](https://img.shields.io/badge/CEM_MoMs-gitee-orange.svg)](https://gitee.com/deltaeecs/CEM_MoMs.git)

## 介绍

提供 CEM_MoMs （主页：[github](https://github.com/deltaeecs/CEM_MoMs.jl), [gitee](https://gitee.com/deltaeecs/CEM_MoMs.jl)） 包的基础信息，包含几何、网格相关类型和函数的定义、接口，基函数相关信息接口。CEM_MoMs 本身被拆分为几个独立的包以方便开发时快速编译， 同时避免在无图形化界面使用时调入绘图相关包而导致报错。

## 安装与测试

### 安装

可直接在 Julia 的包管理模式中安装：

```julia
julia> Pkg.add("MoM_Basics")
```

或

```julia
pkg> add MoM_Basics
```

### 测试

同样可直接在 Julia 的包管理模式中测试包：

```julia
julia> Pkg.test("MoM_Basics")
```

或

```julia
pkg> add MoM_Basics
```
