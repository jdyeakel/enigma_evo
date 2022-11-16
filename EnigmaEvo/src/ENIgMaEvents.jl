abstract type AbstractENIgMaEvent end

struct ColonizationEvent <: AbstractENIgMaEvent
    colonizerId::Int
end

struct PrimaryExtinctionEvent <: AbstractENIgMaEvent
    extId::Int
end

struct SecondaryExtinctionEvent <: AbstractENIgMaEvent
    extId::Int
end

struct ObjectExtinctionEvent <: AbstractENIgMaEvent
    extId::Int
end

struct GlobalExtinctionEvent <: AbstractENIgMaEvent
    extId::Int
end

struct MutationEvent <: AbstractENIgMaEvent
    parentSpecId::Int
    interactentId::Int
    newSpecId::Int
    oldInteraction::InteractionType
    newInteraction::InteractionType
    changeInInteraction::Bool
end

isMutationType(event::AbstractENIgMaEvent, varArgs...) = false

function isMutationType(event::MutationEvent, oldInteraction, newInteraction, changeInInteraction)
    ret = event.oldInteraction == oldInteraction
    ret &= event.newInteraction == newInteraction
    ret &= event.changeInInteraction == changeInInteraction
end

function isMutationType(event::MutationEvent, oldInteraction, newInteraction)
    ret = event.oldInteraction == oldInteraction
    ret &= event.newInteraction == newInteraction 
end