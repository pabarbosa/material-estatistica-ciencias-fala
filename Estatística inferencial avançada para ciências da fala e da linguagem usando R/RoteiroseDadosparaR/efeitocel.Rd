\name{efeitocel}
\alias{efeitocel}
\docType{data}
\title{
Efeito do celular para a medida das frequências dos dois primeiros formantes de vogais orais
}
\description{
Dados de frequências dos dois primeiros formantes de vogais orais em gravação direta e via celular GSM de um informante paulista masculino.
}
\usage{data(efeitocel)}
\format{
  Um data frame com 478 observações de quatro variáveis.
  \describe{
    \item{\code{cond}}{fator com 2 níveis identificando a gravação direta (\code[dir]) e a via celular (\code[cel]).}
    \item{\code{vogal}}{fator com 10 níveis, indicando o monotongo oral em que as frequências de formantes foram medidas.} 
    \item{\code{F1}}{vetor numérico, indicando a frequência do primeiro formante em Hertz da vogal corrrespondente.} 
    \item{\code{F2}}{vetor numérico, indicando a frequência do segundo formante em Hertz da vogal corrrespondente.} 
}
\details{
O informante conversou via celular com a experimentadora e uma placa permitiu gravar seu sinal no celular da experimentadora. Um microfone próximo a sua boca permitiu gravar ao mesmo tempo. sua fala sem intervenção do aparelho. Assim, a única diferença entre os valores de frequências dos formantes entre os sinais de fala é devida ao efeito do filtro do celular.
}
\source{
Dados coletados por Renata Regina Passetti em 2014 no Grupo de Estudos de Prosódia da Fala da Universidade Estadual de Campinas para seu mestrado intitulado "O efeito do telefone celular no sinal da fala: uma análise fonético-acústica com implicações para a verificação de locutor em português brasileiro" e defendido em 2015 na mesma instituição.
}
\references{
}
\examples{
data(efeitocel)
}
\keyword{datasets}
